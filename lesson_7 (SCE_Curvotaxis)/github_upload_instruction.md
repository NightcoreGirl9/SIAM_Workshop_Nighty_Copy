# Uploading Your Final Curvotaxis Project to GitHub

In this final step, you will create your own GitHub repository and upload the completed project files.

The files you need to upload are:

```text
forces.h
SCE_Curvotaxis.cpp
plot_animation_gif.gnu
plot_animation.gnu
```

You may choose to make your repository **public** or **private**.

- **Public** means anyone can view your code.
- **Private** means only you and people you invite can view your code.

For this workshop, either option is fine.

---

## Part 1: Create a New Repository on GitHub

1. Go to [GitHub](https://github.com/) and sign in.

2. In the upper-right corner, click the **+** button.

3. Select **New repository**.

4. Give your repository a name.

   Suggested names:

   ```text
   Curvotaxis_Project
   ```

   or

   ```text
   SCE_Curvotaxis_Model
   ```

5. Choose whether the repository should be **Public** or **Private**.

6. Check the box that says:

   ```text
   Add a README file
   ```

7. Click **Create repository**.

You now have an empty GitHub repository for your final project.

---

## Part 2: Upload the Project Files Using the GitHub Website

This is the easiest method.

1. Open your new repository on GitHub.

2. Click the **Add file** button.

3. Select **Upload files**.

4. Drag and drop the following files into the upload area:

   ```text
   forces.h
   SCE_Curvotaxis.cpp
   plot_animation_gif.gnu
   plot_animation.gnu
   ```

5. Scroll down to the **Commit changes** section.

6. In the commit message box, write something like:

   ```text
   Upload final curvotaxis project files
   ```

7. Click **Commit changes**.

Your files are now uploaded to your GitHub repository.

---

## Part 3: Check That Your Files Uploaded Correctly

After uploading, your repository should contain at least these files:

```text
forces.h
SCE_Curvotaxis.cpp
plot_animation_gif.gnu
plot_animation.gnu
README.md
```

Click each file to make sure GitHub displays the correct contents.

---

## Part 4: Update Your README

Your repository should have a `README.md` file. This file explains what your project does.

At minimum, your README should include:

1. A short project title.
2. A short description of curvotaxis.
3. A short explanation of what the simulation does.
4. Instructions for compiling and running the code.
5. Instructions for plotting the output using Gnuplot.
6. A brief note about what parameters students can modify.

A template will be provided next week.
---

## Optional: Upload Files Using Git Commands

If you already know how to use Git from the terminal, you can also upload the files using the command line.

First, clone your new repository:

```bash
git clone https://github.com/YOUR_USERNAME/YOUR_REPOSITORY_NAME.git
```

Move into the repository folder:

```bash
cd YOUR_REPOSITORY_NAME
```

Copy these files into the folder:

```text
forces.h
SCE_Curvotaxis.cpp
plot_animation_gif.gnu
plot_animation.gnu
```

Then add, commit, and push:

```bash
git add forces.h SCE_Curvotaxis.cpp plot_animation_gif.gnu plot_animation.gnu README.md
git commit -m "Upload final curvotaxis project"
git push
```

Refresh your GitHub page to confirm that the files are there.

---

## Final Checklist

Before submitting your GitHub link, make sure:

- [ ] Your repository exists.
- [ ] Your repository is either public or private, depending on your preference.
- [ ] The file `forces.h` is uploaded.
- [ ] The file `SCE_Curvotaxis.cpp` is uploaded.
- [ ] The file `plot_animation_gif.gnu` is uploaded.
- [ ] The file `plot_animation.gnu` is uploaded.
- [ ] Your `README.md` explains the project.
- [ ] Your code compiles.
- [ ] Your Gnuplot script runs.
- [ ] You copied the repository link to submit or share.

---

## How to Share Your Repository Link

Open your repository on GitHub and copy the URL from your browser.

It should look something like:

```text
https://github.com/YOUR_USERNAME/Curvotaxis_Project
```

Submit or share this link when asked.
