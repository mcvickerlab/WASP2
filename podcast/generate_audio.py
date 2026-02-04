#!/usr/bin/env python3
"""
Generate audio for WASP's Nest podcast episodes.

Supports two TTS backends:
1. ElevenLabs API (premium quality, requires API key)
2. edge-tts (free fallback)

Usage:
    # With ElevenLabs (set ELEVEN_API_KEY environment variable)
    python generate_audio.py --engine elevenlabs

    # With edge-tts (free, default)
    python generate_audio.py --engine edge-tts

    # Regenerate specific episode
    python generate_audio.py --episode 2

    # Verbose output
    python generate_audio.py --verbose
"""

from __future__ import annotations

import argparse
import asyncio
import functools
import logging
import os
import re
import shutil
import sys
import time
from pathlib import Path
from typing import Iterator

# Configure logging
logging.basicConfig(
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# Directories
SCRIPT_DIR = Path(__file__).parent
CHRONICLES_DIR = SCRIPT_DIR / "chronicles"
AUDIO_DIR = SCRIPT_DIR / "audio"

# ElevenLabs has a 5000 character limit per request
ELEVENLABS_CHAR_LIMIT = 5000

# Timeout for TTS operations (5 minutes)
TTS_TIMEOUT_SECONDS = 300


class AudioGenerationError(Exception):
    """Raised when audio generation fails."""
    pass


def clean_markdown(text: str) -> str:
    """Convert markdown to speakable text optimized for TTS."""
    # Remove YAML front matter
    text = re.sub(r'^---.*?---\s*', '', text, flags=re.DOTALL)

    # Remove code blocks
    text = re.sub(r'```[\s\S]*?```', '', text)

    # Remove inline code
    text = re.sub(r'`[^`]+`', '', text)

    # Remove markdown headers but keep text
    text = re.sub(r'^#{1,6}\s*', '', text, flags=re.MULTILINE)

    # Remove bold/italic markers
    text = re.sub(r'\*\*([^*]+)\*\*', r'\1', text)
    text = re.sub(r'\*([^*]+)\*', r'\1', text)

    # Remove links but keep text
    text = re.sub(r'\[([^\]]+)\]\([^)]+\)', r'\1', text)

    # Remove tables
    text = re.sub(r'\|.*\|', '', text)

    # Remove horizontal rules
    text = re.sub(r'^---+$', '', text, flags=re.MULTILINE)

    # Remove episode metadata section
    text = re.sub(r'## Episode Metadata[\s\S]*$', '', text)

    # Remove illumination references
    text = re.sub(r'See:.*?\.md.*', '', text)

    # Clean up whitespace
    text = re.sub(r'\n{3,}', '\n\n', text)
    text = re.sub(r'[ \t]+', ' ', text)

    return text.strip()


def chunk_text(text: str, max_chars: int = ELEVENLABS_CHAR_LIMIT) -> Iterator[str]:
    """
    Split text into chunks that fit within character limits.

    Splits on sentence boundaries to avoid cutting words.
    """
    if len(text) <= max_chars:
        yield text
        return

    # Split on sentence boundaries
    sentences = re.split(r'(?<=[.!?])\s+', text)
    current_chunk = ""

    for sentence in sentences:
        if len(current_chunk) + len(sentence) + 1 <= max_chars:
            current_chunk = f"{current_chunk} {sentence}".strip()
        else:
            if current_chunk:
                yield current_chunk
            # Handle sentences longer than max_chars
            if len(sentence) > max_chars:
                # Split on word boundaries as fallback
                words = sentence.split()
                current_chunk = ""
                for word in words:
                    if len(current_chunk) + len(word) + 1 <= max_chars:
                        current_chunk = f"{current_chunk} {word}".strip()
                    else:
                        if current_chunk:
                            yield current_chunk
                        current_chunk = word
            else:
                current_chunk = sentence

    if current_chunk:
        yield current_chunk


# Exceptions worth retrying (transient network/server issues)
RETRYABLE_EXCEPTIONS = (ConnectionError, TimeoutError, OSError)


def retry_with_backoff(max_retries: int = 3, base_delay: float = 1.0):
    """
    Decorator for retrying functions with exponential backoff.

    Only retries on transient errors (ConnectionError, TimeoutError, OSError).
    Non-retryable errors (ValueError, AuthenticationError, etc.) fail immediately.
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            last_exception = None
            for attempt in range(max_retries):
                try:
                    return func(*args, **kwargs)
                except RETRYABLE_EXCEPTIONS as e:
                    last_exception = e
                    if attempt < max_retries - 1:
                        delay = base_delay * (2 ** attempt)
                        logger.warning(
                            f"Attempt {attempt + 1} failed: {e}. "
                            f"Retrying in {delay:.1f}s..."
                        )
                        time.sleep(delay)
                except Exception:
                    # Non-retryable error - fail immediately
                    raise
            raise last_exception
        return wrapper
    return decorator


async def generate_with_edge_tts(text: str, output_file: Path) -> None:
    """Generate audio using edge-tts (free Microsoft TTS)."""
    try:
        import edge_tts
    except ImportError:
        raise AudioGenerationError(
            "edge-tts not installed. Install with: pip install edge-tts"
        )

    # Voice configuration for Queen Bee character
    voice = "en-US-AriaNeural"  # Try Aria instead of Jenny
    rate = "+0%"  # Normal rate
    pitch = "+0Hz"  # Normal pitch

    try:
        communicate = edge_tts.Communicate(text, voice, rate=rate, pitch=pitch)
        await asyncio.wait_for(
            communicate.save(str(output_file)),
            timeout=TTS_TIMEOUT_SECONDS
        )
    except asyncio.TimeoutError:
        raise AudioGenerationError(
            f"edge-tts timed out after {TTS_TIMEOUT_SECONDS}s"
        )
    except Exception as e:
        raise AudioGenerationError(f"edge-tts failed: {e}")


@retry_with_backoff(max_retries=3)
def generate_with_elevenlabs(text: str, output_file: Path) -> None:
    """Generate audio using ElevenLabs API (premium quality)."""
    try:
        from elevenlabs.client import ElevenLabs
        from elevenlabs import save
    except ImportError:
        raise AudioGenerationError(
            "elevenlabs not installed. Install with: pip install elevenlabs"
        )

    api_key = os.environ.get("ELEVEN_API_KEY")
    if not api_key:
        raise AudioGenerationError(
            "ELEVEN_API_KEY environment variable not set. "
            "Get your API key from https://elevenlabs.io/app/settings/api-keys"
        )

    client = ElevenLabs(api_key=api_key)

    # Use a warm, professional voice for the Queen Bee character
    voice_id = os.environ.get("ELEVEN_VOICE_ID", "21m00Tcm4TlvDq8ikWAM")

    # Handle text chunking for long content
    chunks = list(chunk_text(text, ELEVENLABS_CHAR_LIMIT))

    if len(chunks) == 1:
        # Single chunk - straightforward
        audio = client.text_to_speech.convert(
            voice_id=voice_id,
            text=text,
            model_id="eleven_multilingual_v2",
            output_format="mp3_44100_128",
            voice_settings={
                "stability": 0.5,
                "similarity_boost": 0.75,
            }
        )
        save(audio, str(output_file))
    else:
        # Multiple chunks - generate and concatenate
        logger.info(f"Text split into {len(chunks)} chunks")
        temp_files = []
        try:
            for i, chunk in enumerate(chunks):
                logger.debug(f"Processing chunk {i + 1}/{len(chunks)}")
                temp_file = output_file.with_suffix(f".part{i}.mp3")
                audio = client.text_to_speech.convert(
                    voice_id=voice_id,
                    text=chunk,
                    model_id="eleven_multilingual_v2",
                    output_format="mp3_44100_128",
                    voice_settings={
                        "stability": 0.5,
                        "similarity_boost": 0.75,
                    }
                )
                save(audio, str(temp_file))
                temp_files.append(temp_file)

            # Concatenate using pydub or ffmpeg
            _concatenate_audio_files(temp_files, output_file)
        finally:
            # Cleanup temp files
            for temp_file in temp_files:
                if temp_file.exists():
                    temp_file.unlink()


def _concatenate_audio_files(input_files: list[Path], output_file: Path) -> None:
    """
    Concatenate multiple audio files into one.

    Attempts pydub first (re-encodes at 128kbps), falls back to ffmpeg
    concat filter (stream copy, no re-encoding) if pydub unavailable.
    """
    try:
        from pydub import AudioSegment
        combined = AudioSegment.empty()
        for f in input_files:
            combined += AudioSegment.from_mp3(str(f))
        combined.export(str(output_file), format="mp3", bitrate="128k")
    except ImportError:
        # Fallback to ffmpeg
        import subprocess

        # Try to find ffmpeg
        ffmpeg_cmd = shutil.which("ffmpeg")
        if not ffmpeg_cmd:
            # Try static-ffmpeg package
            try:
                import static_ffmpeg
            except ImportError:
                pass  # Package not installed - acceptable
            else:
                try:
                    static_ffmpeg.add_paths()
                    ffmpeg_cmd = shutil.which("ffmpeg")
                except Exception as e:
                    logger.warning(f"static_ffmpeg.add_paths() failed: {e}")

        if not ffmpeg_cmd:
            raise AudioGenerationError(
                "ffmpeg not found for audio concatenation. Install with:\n"
                "  pip install pydub  (preferred)\n"
                "  pip install static-ffmpeg\n"
                "  conda install -c conda-forge ffmpeg"
            )

        # Create concat file list
        list_file = output_file.with_suffix(".txt")
        with open(list_file, "w") as f:
            for input_file in input_files:
                f.write(f"file '{input_file}'\n")

        try:
            result = subprocess.run(
                [ffmpeg_cmd, "-y", "-f", "concat", "-safe", "0",
                 "-i", str(list_file), "-c", "copy", str(output_file)],
                capture_output=True,
                text=True,
                check=True
            )
        except subprocess.CalledProcessError as e:
            stderr = e.stderr if e.stderr else "Unknown error"
            raise AudioGenerationError(f"Failed to concatenate audio: {stderr}")
        finally:
            if list_file.exists():
                list_file.unlink()


async def generate_episode_audio(
    episode_file: Path,
    output_file: Path,
    engine: str = "edge-tts"
) -> Path:
    """Generate audio for a single episode."""
    logger.info(f"Processing: {episode_file.name}")
    logger.info(f"  Engine: {engine}")

    # Validate input file exists
    if not episode_file.exists():
        raise AudioGenerationError(f"Episode file not found: {episode_file}")

    # Read and clean the markdown
    try:
        content = episode_file.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as e:
        raise AudioGenerationError(f"Failed to read {episode_file}: {e}")

    text = clean_markdown(content)

    # Validate text is not empty
    if not text or len(text.strip()) < 10:
        raise AudioGenerationError(
            f"Episode {episode_file.name} has no speakable content after cleaning"
        )

    logger.debug(f"  Text length: {len(text)} characters")

    # Generate audio based on engine choice
    if engine == "elevenlabs":
        generate_with_elevenlabs(text, output_file)
    else:
        await generate_with_edge_tts(text, output_file)

    # Validate output was created
    if not output_file.exists():
        raise AudioGenerationError(f"Output file was not created: {output_file}")

    file_size = output_file.stat().st_size
    if file_size < 1000:  # Less than 1KB is suspicious
        raise AudioGenerationError(
            f"Output file is too small ({file_size} bytes), generation may have failed"
        )

    logger.info(f"  -> Saved: {output_file.name} ({file_size / 1024:.1f} KB)")
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


async def main() -> int:
    """Main entry point. Returns exit code."""
    parser = argparse.ArgumentParser(
        description="Generate podcast audio from episode scripts"
    )
    parser.add_argument(
        "--engine",
        choices=["edge-tts", "elevenlabs"],
        default="edge-tts",
        help="TTS engine to use (default: edge-tts)"
    )
    parser.add_argument(
        "--episode",
        type=validate_episode_number,
        help="Generate only specific episode number (1-999)"
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

    # Ensure output directory exists
    AUDIO_DIR.mkdir(exist_ok=True)

    # Find episode files
    if args.episode:
        pattern = f"episode-{args.episode:03d}-*.md"
        episodes = list(CHRONICLES_DIR.glob(pattern))
        if not episodes:
            logger.error(f"No episode file found matching: {pattern}")
            return 1
    else:
        episodes = sorted(CHRONICLES_DIR.glob("episode-*.md"))

    if not episodes:
        logger.error("No episode files found in %s", CHRONICLES_DIR)
        return 1

    print(f"Found {len(episodes)} episode(s)")
    print(f"Engine: {args.engine}")
    print("-" * 40)

    errors = []
    for episode_file in episodes:
        output_name = episode_file.stem + ".mp3"
        output_file = AUDIO_DIR / output_name

        try:
            await generate_episode_audio(episode_file, output_file, args.engine)
        except AudioGenerationError as e:
            logger.error(f"Failed to generate {episode_file.name}: {e}")
            errors.append((episode_file.name, str(e)))
        except Exception as e:
            logger.exception(f"Unexpected error processing {episode_file.name}")
            errors.append((episode_file.name, str(e)))

    print("-" * 40)

    if errors:
        print(f"Completed with {len(errors)} error(s):")
        for name, error in errors:
            print(f"  - {name}: {error}")
        return 1

    print("Done! Audio files generated in:", AUDIO_DIR)
    return 0


if __name__ == "__main__":
    sys.exit(asyncio.run(main()))
