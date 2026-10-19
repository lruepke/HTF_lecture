# HTF Lecture

[![website](https://github.com/lruepke/HTF_lecture/actions/workflows/doc.yml/badge.svg)](https://github.com/lruepke/HTF_lecture/actions/workflows/doc.yml)

Lecture notes for the numerical geophysics courses taught in the geophysics master at Kiel University, Germany. Published at **https://lruepke.github.io/HTF_lecture/**.

## Content

Each folder is one course, named after its semester. There are two lecture series:

| Course | Topics | Semesters |
|---|---|---|
| HydrothermalFoam Lecture | Finite volumes, OpenFOAM and HydrothermalFoam, porous media and hydrothermal flow | Winter2020, Summer2022, Winter2023, Summer2025, Dome2025, Winter2026 |
| Finite Element Modeling in Geodynamics | PDEs, finite differences and finite elements in Python, heat diffusion, Stokes flow, folding | Summer2021, Winter2022, Summer2024, Winter2025 |

The newest course of each series is the one to copy from when preparing a new semester.

## Repository layout

```
<Semester>/            one course (Sphinx project)
  Makefile             `make html` builds the course into build/html
  source/conf.py       Sphinx configuration (title, links, logo)
  source/index.rst     table of contents
  source/...           lecture pages (.rst), notebooks, Python code, OpenFOAM cases
home/                  Hugo homepage with one tile per course
  config.toml          homepage settings and the "Lectures" menu
  content/<Semester>/  homepage tile of a course
.github/workflows/doc.yml   builds and deploys the website
copy_lecture.py        creates a new course from an existing one
environment.yaml       conda environment (Python + Poetry)
pyproject.toml         Python packages, managed by Poetry
poetry.lock            exact package versions (used locally and on GitHub)
```

## Python environment

Python itself comes from conda, all packages are managed by Poetry. The same `poetry.lock` is used locally and for the website build, so both use identical package versions.

Set up once:

```bash
conda env create -f environment.yaml   # creates py312_teaching (Python 3.12, Poetry, pandoc)
conda activate py312_teaching
poetry install --no-root               # installs the packages from poetry.lock into this env
```

Add or update a package:

```bash
conda activate py312_teaching
poetry add <package>          # add a new package
poetry update <package>       # update one package within the limits in pyproject.toml
```

Commit `pyproject.toml` and `poetry.lock` afterwards, otherwise the website build will not have the package.

Students do not use this environment. In the FEM courses they install a smaller conda environment, described in `<Semester>/source/general-info/installation.rst`.

## Building the website

### Locally

```bash
conda activate py312_teaching
cd Winter2026
make html                     # result in Winter2026/build/html/index.html
```

For the homepage, run `hugo server` in `home/` and open the address it prints.

Notebooks are not executed during the build, so save them with their outputs.

### On GitHub

Every push to `main` runs `.github/workflows/doc.yml`:

1. Installs Python 3.12, pandoc and the packages from `poetry.lock` (this includes Hugo).
2. Builds the homepage in `home/` with Hugo.
3. Builds every course listed in `lectures="..."`. A course is only rebuilt if its folder changed since its last build (stored in `commit.ID` on the `gh-pages` branch); otherwise the previous pages are reused.
4. Publishes the result to the `gh-pages` branch, which GitHub Pages serves at https://lruepke.github.io/HTF_lecture/. The site can take a few minutes to update after the run.

The build status is shown under the repository's *Actions* tab.

## Creating a new lecture

```bash
python copy_lecture.py <SourceSemester> <NewSemester> "short description for the homepage tile"
# example:
python copy_lecture.py Summer2025 Winter2026 "OpenFOAM, HydrothermalFoam, porous media flow."
```

The script

- copies only files tracked by git (no build output, simulation results or caches),
- sets the "Edit on GitHub" link and copyright year in the new `source/conf.py`,
- adds the course to the homepage menu (`home/config.toml`) and creates its tile (`home/content/<NewSemester>/`),
- adds the course to the `lectures` list in `.github/workflows/doc.yml`.

Then:

1. Update the content by hand: the title (`project`) in `source/conf.py`, dates and course details in `source/general-info/`, and the lecture pages.
2. Test the build with `make html` in the new folder.
3. Commit and push. The new course appears at `https://lruepke.github.io/HTF_lecture/<newsemester>/` (lower case).
