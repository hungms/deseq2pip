## Boot
# devtools::document()
# roxygen2::roxygenise()
# devtools::check()
# devtools::install()
# devtools::build()
# usethis::use_git()
# pkgdown::build_site()
# usethis::use_github_action()
# usethis::use_pkgdown_github_pages()
# usethis::use_bioc_github_action()


## one time
# pak::pak("r-lib/httr2")
# https://stackoverflow.com/questions/72189273/r-cmd-check-github-actions-workflow-failing-on-warnings-notes
# - uses: r-lib/actions/check-r-package@v2
#        with:
#          build_args: 'c("--no-manual", "--no-build-vignettes")'
#          error-on: '"error"'


# Create a branch:
# git checkout -b feature/add-login

# Work and commit incrementally:
# git add .
# git commit -m "feat: Add login form"
# git commit -m "fix: Correct validation logic"
# git commit -m "refactor: Simplify login controller"

# Rebase and squash:
# git rebase -i main
# git rebase -i main

# Force-push the cleaned branch:
# git push origin feature/add-login --force

# Open a pull request and merge:
# Use "Squash and Merge" to combine all changes into one commit on the main branch.

# Tag the release:
# git tag -a v1.0.0 -m "Initial release with login functionality"
# git push origin v1.0.0