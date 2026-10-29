# Publishing TerravAI (package terravai): GitHub, then CRAN

## 1. Put the package on GitHub

```bash
cd terravai
Rscript tools/set_github_user.R YOUR-USERNAME     # fills URL, BugReports, badge, citation
git add -A && git commit -m "Set GitHub URLs"

# Create the repository and push (GitHub CLI):
gh repo create YOUR-USERNAME/terravai --public --source=. --push
# or, without the CLI: create an empty repo named "terravai" at https://github.com/new, then
git remote add origin https://github.com/YOUR-USERNAME/terravai.git
git push -u origin main
```

The included GitHub Actions workflow runs `R CMD check` on Linux, macOS and
Windows on every push. Wait for it to pass (green badge in the README).

## 2. Final checks before CRAN

```r
# install.packages(c("devtools", "urlchecker", "spelling"))
devtools::check(remote = TRUE, manual = TRUE)   # expect 0 errors, 0 warnings
urlchecker::url_check()                         # all URLs must resolve
spelling::spell_check_package()
devtools::check_win_devel()                     # results emailed in ~30 min
devtools::check_win_release()
```

Update `cran-comments.md` with the results.

## 3. Submit to CRAN

```r
devtools::submit_cran()
```

or upload the tarball from `R CMD build .` at the CRAN web form:
https://cran.r-project.org/submit.html

Confirm the email CRAN sends to the maintainer address. New packages are
reviewed by a person; expect replies within a few days to a few weeks.

## 4. After acceptance

* Tag the release on GitHub: `git tag v0.1.0 && git push --tags`
* Add the CRAN badge to the README:
  `[![CRAN status](https://www.r-pkg.org/badges/version/terravai)](https://CRAN.R-project.org/package=terravai)`
* Archive the release on Zenodo for a DOI to cite in papers.
