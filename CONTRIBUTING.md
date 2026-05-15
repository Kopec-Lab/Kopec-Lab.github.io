# Contributing to the Kopec Lab Website

Welcome to the Kopec Lab website repository. This guide explains how to make changes to the site — please read it before pushing anything.

---

## Prerequisites

Make sure you have the following installed on your machine:

- Git
- Ruby (via rbenv — see below)
- Bundler

### First-time setup

```bash
# Install rbenv and Ruby
brew install rbenv ruby-build
rbenv install 3.2.2
rbenv global 3.2.2

# Add rbenv to your shell
echo 'export PATH="$HOME/.rbenv/bin:$PATH"' >> ~/.zshrc
echo 'eval "$(rbenv init -)"' >> ~/.zshrc
source ~/.zshrc

# Install dependencies
gem install bundler
bundle install
```

### Clone the repo

```bash
git clone https://github.com/Kopec-Lab/Kopec-Lab.github.io.git
cd Kopec-Lab.github.io
```

---

## Previewing the site locally

Always preview your changes locally before pushing:

```bash
bundle exec jekyll serve
```

Open `http://localhost:4000` in your browser. The site rebuilds automatically as you save files.

---

## Making changes — always use a branch

The `main` branch is protected and deploys directly to [kopec-lab.com](https://www.kopec-lab.com). **Never push directly to `main`.**

### Workflow

```bash
# 1. Make sure you have the latest version
git checkout main
git pull

# 2. Create a new branch named after what you're changing
git checkout -b your-name/what-you-changed
# e.g. git checkout -b sam/add-profile-page
# e.g. git checkout -b ivan/update-publications

# 3. Make your changes and preview locally

# 4. Stage and commit
git add .
git commit -m "Short description of what you changed"

# 5. Push your branch
git push origin your-name/what-you-changed

# 6. Open a Pull Request on GitHub
# Go to github.com/Kopec-Lab/Kopec-Lab.github.io
# Click "Compare & pull request"
# Add a short description and request review from @wojtekkopec
```

Wojciech will review and merge. The site redeploys automatically once merged — allow ~5 minutes.

---

## What lives where

| What | Where |
|---|---|
| Front page content | `_pages/about.md` |
| People page | `_pages/profiles.md` + `_pages/about_<name>.md` |
| Profile images | `assets/img/` |
| Blog posts | `_posts/YYYY-MM-DD-title.md` |
| News announcements | `_news/YYYY-MM-DD-title.md` |
| Publications | `_bibliography/papers.bib` |
| Research projects | `_projects/` |
| Teaching content | `_pages/teaching.md` |
| Gallery images | `assets/img/gallery/` + `_data/gallery.yml` |
| Site configuration | `_config.yml` |

---

## Adding a blog post

1. Create a new file in `_posts/` named `YYYY-MM-DD-short-title.md`
2. Use this frontmatter template:

```yaml
---
layout: post
title: Your Post Title
date: YYYY-MM-DD 12:00:00+0000
description: One sentence summary shown in the post list.
tags: [tag1, tag2]
categories: science
related_posts: false
toc:
  sidebar: left
---
```

3. Write your content below the `---` in Markdown
4. Images go in `assets/img/posts/` and are inserted with:

```liquid
{% include figure.liquid path="assets/img/posts/your-image.png" caption="Your caption." zoomable=true %}
```

---

## Adding a news announcement

Create a file in `_news/` named `YYYY-MM-DD-short-title.md`:

```yaml
---
layout: post
date: YYYY-MM-DD 00:00:00+0000
inline: true
related_posts: false
---

Your announcement text here. Keep it to 1-3 sentences.
```

---

## Adding yourself to the People page

1. Add your photo to `assets/img/` (JPEG or PNG, ideally square, under 1 MB)
2. Create `_pages/about_<yourname>.md` with a short bio in Markdown
3. Add an entry to `_pages/profiles.md` under the `profiles:` list:

```yaml
  - align: left
    image: yourname.jpg
    content: about_yourname.md
    image_circular: false
```

4. Open a pull request as usual

---

## Things to avoid

- Do not commit large binary files (trajectories, raw simulation data, videos over 5 MB)
- Do not commit `.DS_Store` or other system files (already in `.gitignore`)
- Do not edit `_config.yml` unless you know what you're changing — it affects the whole site
- Do not push directly to `main`

---

## Questions?

Ask Wojciech or open a GitHub Issue in this repository.
