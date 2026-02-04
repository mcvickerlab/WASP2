# Illuminations - Visual Diagrams

This directory contains Mermaid diagrams and visual aids for podcast episodes.

## Purpose

Illuminations are visual companions to Buzz Reports, helping illustrate:
- Architecture changes
- Feature workflows
- Data flow diagrams
- Version comparisons

## File Naming Convention

```
illumination-{episode_number}-{topic}.md
```

Example: `illumination-001-new-counting-module.md`

## Template

```markdown
# Illumination: [Topic]
# Episode: [NUMBER]

## Overview Diagram

\`\`\`mermaid
graph TD
    A[Input] --> B[Process]
    B --> C[Output]
\`\`\`

## Detailed Flow

\`\`\`mermaid
sequenceDiagram
    participant User
    participant WASP2
    participant Output
    User->>WASP2: Run analysis
    WASP2->>Output: Generate results
\`\`\`
```

## Rendering

Diagrams can be rendered using:
- Mermaid CLI: `mmdc -i input.md -o output.png`
- GitHub's built-in Mermaid support
- VS Code Mermaid extensions
