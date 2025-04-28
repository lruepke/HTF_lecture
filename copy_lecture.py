import shutil
import os
import re
import sys
from datetime import datetime

def copy_lecture(source_dir, destination_dir):
    if not os.path.exists(source_dir):
        print(f"Source directory '{source_dir}' does not exist.")
        return False
    if os.path.exists(destination_dir):
        print(f"Destination directory '{destination_dir}' already exists.")
        return False

    shutil.copytree(source_dir, destination_dir)
    print(f"Lecture copied from '{source_dir}' to '{destination_dir}' successfully.")
    return True

def update_config(config_path, new_name):
    with open(config_path, 'r') as file:
        lines = file.readlines()

    lecture_exists = any(f'name = "{new_name.replace("_", " ")}"' in line or f'name = "{new_name}"' in line for line in lines)
    
    if lecture_exists:
        print(f"Lecture '{new_name}' already exists in config.toml. Skipping adding it again.")
        return

    last_weight = 0
    for line in lines:
        if "weight" in line and "Lectures" in "".join(lines[max(0, lines.index(line)-5):lines.index(line)]):
            try:
                last_weight = max(last_weight, int(line.strip().split("=")[-1]))
            except:
                pass

    new_weight = last_weight + 1

    insert_index = None
    for i in range(len(lines)-1, -1, -1):
        if 'parent = "Lectures"' in lines[i]:
            insert_index = i
            break

    if insert_index is None:
        print("Could not find 'Lectures' section in config.toml!")
        return

    new_block = f"""
[[Languages.en.menu.main]]
parent = "Lectures"
name = "{new_name.replace('_', ' ')}"
url = "{new_name.lower()}"
weight = {new_weight}
"""

    lines.insert(insert_index + 5, new_block)

    with open(config_path, 'w') as file:
        file.writelines(lines)

    print(f"Added new lecture '{new_name}' to {config_path} with weight {new_weight}.")

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

def create_homepage_tile(new_name):
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
description: "Lecture notes for {pretty_title}: Advanced OpenFOAM, fluid modeling, geophysical flows."
type: "docs"
color: "73, 147, 65"
---

<script type="text/javascript">
  window.open("/{new_name.lower()}/", "_self");
</script>
"""

    with open(index_md_path, 'w') as f:
        f.write(index_md_content)

    print(f"Created homepage tile at {index_md_path}.")

if __name__ == "__main__":
    default_source = "Winter2023"
    default_destination = "Summer2025"

    if len(sys.argv) == 3:
        source = sys.argv[1]
        destination = sys.argv[2]
        print(f"Using provided arguments: {source} -> {destination}")
    else:
        source = default_source
        destination = default_destination
        print(f"No arguments provided, using defaults: {source} -> {destination}")

    config_file = os.path.join("home", "config.toml")
    workflow_file = os.path.join(".github", "workflows", "doc.yml")

    if copy_lecture(source, destination):
        update_config(config_file, destination)
        update_workflow(workflow_file, destination)
        create_homepage_tile(destination)
