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
"""

import argparse
import asyncio
import os
import re
from pathlib import Path

# Directories
SCRIPT_DIR = Path(__file__).parent
CHRONICLES_DIR = SCRIPT_DIR / "chronicles"
AUDIO_DIR = SCRIPT_DIR / "audio"


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


async def generate_with_edge_tts(text: str, output_file: Path):
    """Generate audio using edge-tts (free Microsoft TTS)."""
    import edge_tts

    # Voice configuration for Queen Bee character
    voice = "en-US-JennyNeural"
    rate = "-5%"  # Slightly slower for clarity
    pitch = "+2Hz"  # Slightly higher for Queen Bee character

    communicate = edge_tts.Communicate(text, voice, rate=rate, pitch=pitch)
    await communicate.save(str(output_file))


def generate_with_elevenlabs(text: str, output_file: Path):
    """Generate audio using ElevenLabs API (premium quality)."""
    from elevenlabs.client import ElevenLabs
    from elevenlabs import save

    api_key = os.environ.get("ELEVEN_API_KEY")
    if not api_key:
        raise ValueError(
            "ELEVEN_API_KEY environment variable not set. "
            "Get your API key from https://elevenlabs.io/app/settings/api-keys"
        )

    client = ElevenLabs(api_key=api_key)

    # Use a warm, professional voice for the Queen Bee character
    # Rachel (21m00Tcm4TlvDq8ikWAM) - warm, professional narrator voice
    # Or use your own cloned voice ID
    voice_id = os.environ.get("ELEVEN_VOICE_ID", "21m00Tcm4TlvDq8ikWAM")

    audio = client.text_to_speech.convert(
        voice_id=voice_id,
        text=text,
        model_id="eleven_multilingual_v2",
        output_format="mp3_44100_128",
        voice_settings={
            "stability": 0.5,  # More expressive
            "similarity_boost": 0.75,  # Natural voice
        }
    )

    save(audio, str(output_file))


async def generate_episode_audio(
    episode_file: Path,
    output_file: Path,
    engine: str = "edge-tts"
):
    """Generate audio for a single episode."""
    print(f"Processing: {episode_file.name}")
    print(f"  Engine: {engine}")

    # Read and clean the markdown
    content = episode_file.read_text()
    text = clean_markdown(content)

    # Generate audio based on engine choice
    if engine == "elevenlabs":
        generate_with_elevenlabs(text, output_file)
    else:
        await generate_with_edge_tts(text, output_file)

    print(f"  -> Saved: {output_file.name}")
    return output_file


async def main():
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
        type=int,
        help="Generate only specific episode number (1, 2, or 3)"
    )
    args = parser.parse_args()

    AUDIO_DIR.mkdir(exist_ok=True)

    # Find episode files
    if args.episode:
        pattern = f"episode-{args.episode:03d}-*.md"
        episodes = list(CHRONICLES_DIR.glob(pattern))
        if not episodes:
            print(f"No episode file found matching: {pattern}")
            return
    else:
        episodes = sorted(CHRONICLES_DIR.glob("episode-*.md"))

    if not episodes:
        print("No episode files found!")
        return

    print(f"Found {len(episodes)} episode(s)")
    print(f"Engine: {args.engine}")
    print("-" * 40)

    for episode_file in episodes:
        # Create output filename
        output_name = episode_file.stem + ".mp3"
        output_file = AUDIO_DIR / output_name

        await generate_episode_audio(episode_file, output_file, args.engine)

    print("-" * 40)
    print("Done! Audio files generated in:", AUDIO_DIR)


if __name__ == "__main__":
    asyncio.run(main())
