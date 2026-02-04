#!/usr/bin/env python3
"""
Audio post-processing pipeline for WASP's Nest podcast.

Applies professional audio enhancement using ffmpeg:
1. Noise reduction (afftdn filter)
2. High-pass filter (remove rumble < 80Hz)
3. Low-pass filter (remove hiss > 12kHz)
4. Compression (reduce dynamic range)
5. Loudness normalization (podcast standard: -16 LUFS)

Requirements:
    - ffmpeg with libavfilter (auto-detects static-ffmpeg if installed)

Usage:
    python enhance_audio.py                    # Enhance all episodes
    python enhance_audio.py --episode 2        # Enhance specific episode
    python enhance_audio.py --dry-run          # Show commands without running
    python enhance_audio.py --verbose          # Verbose output
"""

from __future__ import annotations

import argparse
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator

# Configure logging
logging.basicConfig(
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

SCRIPT_DIR = Path(__file__).parent
AUDIO_DIR = SCRIPT_DIR / "audio"
ENHANCED_DIR = SCRIPT_DIR / "audio_enhanced"

# Processing timeout (10 minutes per file)
PROCESS_TIMEOUT_SECONDS = 600


class AudioEnhanceError(Exception):
    """Raised when audio enhancement fails."""
    pass


def find_ffmpeg() -> str:
    """
    Find ffmpeg executable, trying multiple sources.

    Returns:
        Path to ffmpeg executable

    Raises:
        AudioEnhanceError: If ffmpeg is not found
    """
    # Try system ffmpeg first
    ffmpeg_path = shutil.which("ffmpeg")
    if ffmpeg_path:
        logger.debug(f"Found system ffmpeg: {ffmpeg_path}")
        return ffmpeg_path

    # Try static-ffmpeg package
    try:
        import static_ffmpeg
    except ImportError:
        pass  # Package not installed - try other options
    else:
        # Package is installed - failure here is an error, not a fallback
        try:
            ffmpeg_path, _ = static_ffmpeg.run.get_or_fetch_platform_executables_else_raise()
            logger.debug(f"Found static-ffmpeg: {ffmpeg_path}")
            return ffmpeg_path
        except Exception as e:
            raise AudioEnhanceError(
                f"static-ffmpeg is installed but failed: {e}\n"
                "Try: pip uninstall static-ffmpeg && pip install static-ffmpeg"
            )

    # Try common installation paths
    common_paths = [
        "/usr/bin/ffmpeg",
        "/usr/local/bin/ffmpeg",
        os.path.expanduser("~/.local/bin/ffmpeg"),
    ]
    for path in common_paths:
        if os.path.isfile(path) and os.access(path, os.X_OK):
            logger.debug(f"Found ffmpeg at: {path}")
            return path

    raise AudioEnhanceError(
        "ffmpeg not found. Install with:\n"
        "  pip install static-ffmpeg && python -c 'import static_ffmpeg; static_ffmpeg.add_paths()'\n"
        "  or: conda install -c conda-forge ffmpeg\n"
        "  or: apt-get install ffmpeg"
    )


def build_ffmpeg_filter() -> str:
    """
    Build the ffmpeg audio filter chain for podcast enhancement.

    Filter chain:
    1. afftdn - FFT-based noise reduction (reduces steady background noise)
    2. highpass - Remove low-frequency rumble (< 80Hz)
    3. lowpass - Remove high-frequency hiss (> 12kHz)
    4. acompressor - Dynamic range compression (voice clarity)
    5. loudnorm - EBU R128 loudness normalization (-16 LUFS for podcasts)
    """
    filters = [
        # Noise reduction: removes steady background noise
        # nr=12 = noise reduction in dB, nf=-25 = noise floor
        "afftdn=nr=12:nf=-25",

        # High-pass filter: remove rumble below 80Hz
        # Human voice is 85Hz-255Hz fundamental, so 80Hz is safe
        "highpass=f=80",

        # Low-pass filter: remove hiss above 12kHz
        # Preserves voice clarity while removing high-freq artifacts
        "lowpass=f=12000",

        # Dynamic range compression for consistent volume
        # threshold=-20dB, ratio=4:1, attack=5ms, release=50ms
        "acompressor=threshold=-20dB:ratio=4:attack=5:release=50",

        # Loudness normalization to podcast standard
        # -16 LUFS is the standard for podcasts (Spotify, Apple Podcasts)
        # TP=-1.5 = true peak limit to prevent clipping
        "loudnorm=I=-16:TP=-1.5:LRA=11",
    ]

    return ",".join(filters)


@contextmanager
def temp_file_context(suffix: str = ".mp3") -> Iterator[Path]:
    """Context manager for temporary file with guaranteed cleanup."""
    fd, path = tempfile.mkstemp(suffix=suffix)
    os.close(fd)
    temp_path = Path(path)
    try:
        yield temp_path
    finally:
        if temp_path.exists():
            try:
                temp_path.unlink()
            except OSError as e:
                logger.warning(f"Failed to cleanup temp file {temp_path}: {e}")


def validate_audio_file(path: Path) -> None:
    """
    Validate that a file is a valid audio file.

    Raises:
        AudioEnhanceError: If file is invalid
    """
    if not path.exists():
        raise AudioEnhanceError(f"File not found: {path}")

    if not path.is_file():
        raise AudioEnhanceError(f"Not a file: {path}")

    # Check file size (minimum 1KB for valid audio)
    size = path.stat().st_size
    if size < 1024:
        raise AudioEnhanceError(
            f"File too small ({size} bytes), may be corrupted: {path}"
        )

    # Check file extension
    if path.suffix.lower() not in {".mp3", ".wav", ".m4a", ".ogg", ".flac"}:
        logger.warning(f"Unusual audio extension: {path.suffix}")


def enhance_audio(
    input_file: Path,
    output_file: Path,
    ffmpeg_path: str,
    dry_run: bool = False
) -> Path:
    """
    Apply audio enhancement to a single file.

    Args:
        input_file: Path to input audio file
        output_file: Path for enhanced output
        ffmpeg_path: Path to ffmpeg executable
        dry_run: If True, print command without executing

    Returns:
        Path to the enhanced audio file

    Raises:
        AudioEnhanceError: If enhancement fails
    """
    # Validate input
    validate_audio_file(input_file)

    filter_chain = build_ffmpeg_filter()

    cmd = [
        ffmpeg_path,
        "-y",  # Overwrite output
        "-i", str(input_file),
        "-af", filter_chain,
        "-c:a", "libmp3lame",  # MP3 output
        "-b:a", "192k",  # 192kbps bitrate
        "-ar", "44100",  # 44.1kHz sample rate
        str(output_file)
    ]

    logger.info(f"Processing: {input_file.name}")

    if dry_run:
        print(f"  Command: {' '.join(cmd)}")
        return output_file

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            timeout=PROCESS_TIMEOUT_SECONDS
        )
        logger.debug(f"ffmpeg stdout: {result.stdout}")
    except subprocess.TimeoutExpired:
        raise AudioEnhanceError(
            f"ffmpeg timed out after {PROCESS_TIMEOUT_SECONDS}s for {input_file.name}"
        )
    except subprocess.CalledProcessError as e:
        raise AudioEnhanceError(
            f"ffmpeg failed for {input_file.name}: {e.stderr}"
        )

    # Validate output was created and is valid
    if not output_file.exists():
        raise AudioEnhanceError(f"Output file was not created: {output_file}")

    output_size = output_file.stat().st_size
    input_size = input_file.stat().st_size

    # Output should be reasonably sized (at least 10% of input)
    if output_size < input_size * 0.1:
        raise AudioEnhanceError(
            f"Output file suspiciously small ({output_size} bytes vs "
            f"{input_size} bytes input), enhancement may have failed"
        )

    logger.info(f"  -> Enhanced: {output_file.name} ({output_size / 1024 / 1024:.1f} MB)")
    return output_file


def validate_episode_number(value: str) -> int:
    """Validate episode number is a positive integer."""
    try:
        episode = int(value)
        if episode < 1 or episode > 999:
            raise argparse.ArgumentTypeError(
                f"Episode number must be between 1 and 999, got {episode}"
            )
        return episode
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"Episode must be a number, got '{value}'"
        )


def main() -> int:
    """Main entry point. Returns exit code."""
    parser = argparse.ArgumentParser(
        description="Enhance podcast audio with noise reduction and normalization"
    )
    parser.add_argument(
        "--episode",
        type=validate_episode_number,
        help="Enhance only specific episode number (1-999)"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show ffmpeg commands without running"
    )
    parser.add_argument(
        "--in-place",
        action="store_true",
        help="Overwrite original files instead of creating enhanced copies"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose output"
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug output"
    )
    args = parser.parse_args()

    # Configure logging level
    if args.debug:
        logger.setLevel(logging.DEBUG)
    elif args.verbose:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)

    # Find ffmpeg
    try:
        ffmpeg_path = find_ffmpeg()
    except AudioEnhanceError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    # Determine output directory
    if args.in_place:
        output_dir = AUDIO_DIR
    else:
        output_dir = ENHANCED_DIR
        output_dir.mkdir(exist_ok=True)

    # Find audio files
    if args.episode:
        pattern = f"episode-{args.episode:03d}-*.mp3"
        audio_files = list(AUDIO_DIR.glob(pattern))
        if not audio_files:
            print(f"No audio file found matching: {pattern}", file=sys.stderr)
            return 1
    else:
        audio_files = sorted(AUDIO_DIR.glob("episode-*.mp3"))

    if not audio_files:
        print(f"No audio files found in {AUDIO_DIR}", file=sys.stderr)
        return 1

    print(f"Found {len(audio_files)} audio file(s)")
    print(f"Output: {output_dir}")
    print(f"ffmpeg: {ffmpeg_path}")
    print("-" * 40)

    errors = []
    for audio_file in audio_files:
        try:
            if args.in_place:
                # Create temp file, enhance, then replace original
                # Use manual temp file management to preserve on move failure
                fd, temp_path_str = tempfile.mkstemp(suffix=".mp3")
                os.close(fd)
                temp_file = Path(temp_path_str)
                try:
                    enhance_audio(audio_file, temp_file, ffmpeg_path, args.dry_run)
                    if not args.dry_run:
                        try:
                            shutil.move(str(temp_file), str(audio_file))
                        except Exception as e:
                            # Keep enhanced file for recovery
                            backup_path = audio_file.with_suffix(".enhanced.mp3")
                            shutil.copy(str(temp_file), str(backup_path))
                            raise AudioEnhanceError(
                                f"Failed to replace original: {e}. "
                                f"Enhanced version saved to: {backup_path}"
                            )
                finally:
                    # Only cleanup if file still exists (wasn't moved)
                    if temp_file.exists():
                        try:
                            temp_file.unlink()
                        except OSError:
                            pass
            else:
                output_file = output_dir / audio_file.name
                enhance_audio(audio_file, output_file, ffmpeg_path, args.dry_run)
        except AudioEnhanceError as e:
            logger.error(str(e))
            errors.append((audio_file.name, str(e)))
        except Exception as e:
            logger.exception(f"Unexpected error processing {audio_file.name}")
            errors.append((audio_file.name, str(e)))

    print("-" * 40)

    if errors:
        print(f"Completed with {len(errors)} error(s):")
        for name, error in errors:
            print(f"  - {name}: {error}")
        return 1

    print("Done! Enhanced audio files in:", output_dir)
    print()
    print("Enhancement applied:")
    print("  - Noise reduction (afftdn)")
    print("  - High-pass filter (80Hz)")
    print("  - Low-pass filter (12kHz)")
    print("  - Dynamic compression (4:1 ratio)")
    print("  - Loudness normalization (-16 LUFS)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
