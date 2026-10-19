"""Create a new lecture from an existing one.

Usage:
    python copy_lecture.py <SourceLecture> <NewLecture> ["tile description"]

Example:
    python copy_lecture.py Summer2025 Winter2026 "OpenFOAM, hydrothermal flow, porous media."

Only files tracked by git (plus new files not covered by .gitignore) are copied,
so build output, simulation results, caches and editor folders stay behind.
"""
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime


def list_source_files(source_dir):
    """Files in source_dir that git tracks or would track (respects .gitignore)."""
    result = subprocess.run(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z", "--", source_dir],
        capture_output=True, text=True, check=True,
    )
    files = [f for f in result.stdout.split("\0") if f]
    # git ls-files --cached also lists files that were deleted but not yet committed
    return [f for f in files if os.path.isfile(f)]


def copy_lecture(source_dir, destination_dir):
    if not os.path.isdir(source_dir):
        print(f"Source directory '{source_dir}' does not exist.")
        return False
    if os.path.exists(destination_dir):
        print(f"Destination directory '{destination_dir}' already exists.")
        return False

    files = list_source_files(source_dir)
    for src in files:
        dst = os.path.join(destination_dir, os.path.relpath(src, source_dir))
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copy2(src, dst)

    print(f"Copied {len(files)} files from '{source_dir}' to '{destination_dir}'.")
    return True


def update_sphinx_conf(conf_path, new_name):
    """Point the 'Edit on GitHub' link at the new folder and update the copyright year."""
    with open(conf_path, 'r') as file:
        content = file.read()

    year = datetime.today().year
    replacements = [
        (r"'github_version':\s*'[^']*'", f"'github_version': 'main/{new_name}/'"),
        (r"'conf_py_path':\s*'[^']*'", "'conf_py_path': 'source/'"),
        (r"copyright\s*=\s*'\d{4}", f"copyright = '{year}"),
    ]
    for pattern, replacement in replacements:
        content, count = re.subn(pattern, replacement, content)
        if count == 0:
            print(f"Warning: '{pattern}' not found in {conf_path}, please update it by hand.")

    with open(conf_path, 'w') as file:
        file.write(content)

    print(f"Updated GitHub link and copyright year in {conf_path}.")


def update_config(config_path, new_name):
    """Add the lecture to the 'Lectures' dropdown menu of the Hugo homepage."""
    with open(config_path, 'r') as file:
        lines = file.readlines()

    if any(re.match(rf'\s*url\s*=\s*"{new_name.lower()}"', line) for line in lines):
        print(f"Lecture '{new_name}' already exists in {config_path}. Skipping.")
        return

    # Find all active (not commented out) menu entries with parent = "Lectures"
    last_weight = 0
    insert_index = None
    indent = "  "
    for i, line in enumerate(lines):
        if line.strip() != 'parent = "Lectures"':
            continue
        indent = line[:len(line) - len(line.lstrip())]
        for j in range(i + 1, len(lines)):
            stripped = lines[j].strip()
            if stripped.startswith("[") or stripped.startswith("#"):
                break
            match = re.match(r'weight\s*=\s*(\d+)', stripped)
            if match:
                last_weight = max(last_weight, int(match.group(1)))
                insert_index = j + 1
                break

    if insert_index is None:
        print(f"Could not find the 'Lectures' menu in {config_path}!")
        return

    new_weight = last_weight + 1
    new_block = (
        f"\n{indent}[[Languages.en.menu.main]]\n"
        f'{indent}parent = "Lectures"\n'
        f'{indent}name = "{new_name}"\n'
        f'{indent}url = "{new_name.lower()}"\n'
        f"{indent}weight = {new_weight}\n"
    )
    lines.insert(insert_index, new_block)

    with open(config_path, 'w') as file:
        file.writelines(lines)

    print(f"Added '{new_name}' to {config_path} with weight {new_weight}.")


def update_workflow(workflow_path, new_lecture):
    with open(workflow_path, 'r') as file:
        content = file.read()

    lectures_pattern = r'lectures="([^"]+)"'
    match = re.search(lectures_pattern, content)
    if not match:
        print("Could not find lectures list in workflow file!")
        return

    lectures = match.group(1).split()
    if new_lecture in lectures:
        print(f"{new_lecture} is already in the workflow.")
        return

    lectures.append(new_lecture)
    new_lectures = " ".join(lectures)
    new_content = re.sub(lectures_pattern, f'lectures="{new_lectures}"', content)

    with open(workflow_path, 'w') as file:
        file.write(new_content)

    print(f"Added '{new_lecture}' to the lectures list in {workflow_path}.")


def create_homepage_tile(new_name, description):
    content_dir = os.path.join("home", "content", new_name)
    os.makedirs(content_dir, exist_ok=True)

    index_md_path = os.path.join(content_dir, "_index.en.md")
    if os.path.exists(index_md_path):
        print(f"_index.en.md for '{new_name}' already exists. Skipping creation.")
        return

    pretty_title = re.sub(r'(\D+)(\d+)', r'\1 \2', new_name)

    today = datetime.today().strftime('%Y-%m-%dT00:00:00+00:00')

    index_md_content = f"""---
title: "{pretty_title}"
date: {today}
icon: "fas fa-landmark"
description: "Lecture notes for {pretty_title}: {description}"
type: "docs"
color: "73, 147, 65"
---

<script type="text/javascript">
  window.open("/HTF_lecture/{new_name.lower()}/", "_self");
</script>
"""

    with open(index_md_path, 'w') as f:
        f.write(index_md_content)

    print(f"Created homepage tile at {index_md_path}.")


if __name__ == "__main__":
    if len(sys.argv) not in (3, 4):
        print(__doc__)
        sys.exit(1)

    source = sys.argv[1].rstrip("/")
    destination = sys.argv[2].rstrip("/")
    description = sys.argv[3] if len(sys.argv) == 4 else "TODO: describe the lecture."

    # Paths below are relative to the repository root
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    config_file = os.path.join("home", "config.toml")
    workflow_file = os.path.join(".github", "workflows", "doc.yml")
    conf_file = os.path.join(destination, "source", "conf.py")

    if copy_lecture(source, destination):
        update_sphinx_conf(conf_file, destination)
        update_config(config_file, destination)
        update_workflow(workflow_file, destination)
        create_homepage_tile(destination, description)
