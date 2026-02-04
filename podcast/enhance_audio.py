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
    - ffmpeg with libavfilter
    - pydub (optional, for format conversion)

Usage:
    python enhance_audio.py                    # Enhance all episodes
    python enhance_audio.py --episode 2        # Enhance specific episode
    python enhance_audio.py --dry-run          # Show commands without running
"""

import argparse
import subprocess
import shutil
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
AUDIO_DIR = SCRIPT_DIR / "audio"
ENHANCED_DIR = SCRIPT_DIR / "audio_enhanced"


def check_ffmpeg():
    """Verify ffmpeg is available."""
    if not shutil.which("ffmpeg"):
        raise RuntimeError(
            "ffmpeg not found. Install with:\n"
            "  conda install -c conda-forge ffmpeg\n"
            "  or: apt-get install ffmpeg"
        )


def build_ffmpeg_filter():
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


def enhance_audio(
    input_file: Path,
    output_file: Path,
    dry_run: bool = False
) -> Path:
    """
    Apply audio enhancement to a single file.

    Args:
        input_file: Path to input audio file
        output_file: Path for enhanced output
        dry_run: If True, print command without executing

    Returns:
        Path to the enhanced audio file
    """
    filter_chain = build_ffmpeg_filter()

    cmd = [
        "ffmpeg",
        "-y",  # Overwrite output
        "-i", str(input_file),
        "-af", filter_chain,
        "-c:a", "libmp3lame",  # MP3 output
        "-b:a", "192k",  # 192kbps bitrate
        "-ar", "44100",  # 44.1kHz sample rate
        str(output_file)
    ]

    print(f"Processing: {input_file.name}")

    if dry_run:
        print(f"  Command: {' '.join(cmd)}")
        return output_file

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        print(f"  -> Enhanced: {output_file.name}")
        return output_file
    except subprocess.CalledProcessError as e:
        print(f"  Error: {e.stderr}")
        raise


def main():
    parser = argparse.ArgumentParser(
        description="Enhance podcast audio with noise reduction and normalization"
    )
    parser.add_argument(
        "--episode",
        type=int,
        help="Enhance only specific episode number (1, 2, or 3)"
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
    args = parser.parse_args()

    check_ffmpeg()

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
            print(f"No audio file found matching: {pattern}")
            return
    else:
        audio_files = sorted(AUDIO_DIR.glob("episode-*.mp3"))

    if not audio_files:
        print("No audio files found in", AUDIO_DIR)
        return

    print(f"Found {len(audio_files)} audio file(s)")
    print(f"Output: {output_dir}")
    print("-" * 40)

    for audio_file in audio_files:
        if args.in_place:
            # Create temp file, then replace original
            temp_file = audio_file.with_suffix(".enhanced.mp3")
            enhance_audio(audio_file, temp_file, args.dry_run)
            if not args.dry_run:
                temp_file.replace(audio_file)
        else:
            output_file = output_dir / audio_file.name
            enhance_audio(audio_file, output_file, args.dry_run)

    print("-" * 40)
    print("Done! Enhanced audio files in:", output_dir)
    print()
    print("Enhancement applied:")
    print("  - Noise reduction (afftdn)")
    print("  - High-pass filter (80Hz)")
    print("  - Low-pass filter (12kHz)")
    print("  - Dynamic compression (4:1 ratio)")
    print("  - Loudness normalization (-16 LUFS)")


if __name__ == "__main__":
    main()
