#!/usr/bin/env python3
"""Generate audio for WASP's Nest podcast episodes using edge-tts."""

import asyncio
import re
from pathlib import Path
import edge_tts

# Voice configuration
VOICE = "en-US-JennyNeural"
RATE = "-5%"  # Slightly slower for clarity
PITCH = "+2Hz"  # Slightly higher for Queen Bee character

# Directories
SCRIPT_DIR = Path(__file__).parent
CHRONICLES_DIR = SCRIPT_DIR / "chronicles"
AUDIO_DIR = SCRIPT_DIR / "audio"


def clean_markdown(text: str) -> str:
    """Convert markdown to speakable text."""
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

    # Convert emotion tags to pauses
    text = re.sub(r'\[happy buzz\]', '...', text)
    text = re.sub(r'\[excited waggle\]', '...', text)
    text = re.sub(r'\[pause\]', '...', text)
    text = re.sub(r'\[.*?tone\]', '', text)
    text = re.sub(r'\[.*?buzz\]', '...', text)
    text = re.sub(r'\[.*?\]', '', text)

    # Remove episode metadata section
    text = re.sub(r'## Episode Metadata[\s\S]*$', '', text)

    # Remove illumination references
    text = re.sub(r'See:.*?\.md.*', '', text)

    # Clean up whitespace
    text = re.sub(r'\n{3,}', '\n\n', text)
    text = re.sub(r'[ \t]+', ' ', text)

    return text.strip()


async def generate_episode_audio(episode_file: Path, output_file: Path):
    """Generate audio for a single episode."""
    print(f"Processing: {episode_file.name}")

    # Read and clean the markdown
    content = episode_file.read_text()
    text = clean_markdown(content)

    # Generate audio
    communicate = edge_tts.Communicate(text, VOICE, rate=RATE, pitch=PITCH)
    await communicate.save(str(output_file))

    print(f"  -> Saved: {output_file.name}")


async def main():
    """Generate audio for all episodes."""
    AUDIO_DIR.mkdir(exist_ok=True)

    # Find episode files
    episodes = sorted(CHRONICLES_DIR.glob("episode-*.md"))

    if not episodes:
        print("No episode files found!")
        return

    print(f"Found {len(episodes)} episodes")
    print(f"Voice: {VOICE}")
    print("-" * 40)

    for episode_file in episodes:
        # Create output filename
        output_name = episode_file.stem.replace("episode-", "episode-") + ".mp3"
        output_file = AUDIO_DIR / output_name

        await generate_episode_audio(episode_file, output_file)

    print("-" * 40)
    print("Done! Audio files generated in:", AUDIO_DIR)


if __name__ == "__main__":
    asyncio.run(main())
