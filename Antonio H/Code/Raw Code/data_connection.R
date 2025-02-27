### Connecting to GitHub, not sure but hopefully it work
## Antonio H.
# 2/18/25

# the packages I know I need for sure
install.packages("usethis")
install.packages("gitcreds")

#
library(gitcreds)
library(usethis)

# connecting to my system
usethis::use_git_config(user.name = "ADHolmes1999", user.email = "antonioholmes0@gmail.com")

# creating a token
usethis::create_github_token()

# getting connected
gitcreds::gitcreds_set()

# I ran into a problem so I will download Git on my computer then connect to it using Tools >Global options > Git/SVN tab

# verifying it worked this time
usethis::git_sitrep()

# cloning my Git Repo
usethis::create_from_github("https://github.com/ADHolmes1999/LBSL.git")

# initializing git (well... that's what I think it does)
usethis::use_git()

# all best to check status
system("git status")

# always pull before pushing
system("git pull origin main")

# check status one mo'gen
system("git status")

# If I ever need to push, here is all the code
system("git add .")
system("git commit")
system("git push orign main")