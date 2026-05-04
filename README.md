Hi all,

You can install the **NetClusterSim** R package directly from our private GitHub repository. Please follow the steps below.

### 1. Install required package

If you don’t already have it:

```r
install.packages("remotes")
```

### 2. Set up GitHub access (one-time)

Since the repository is private, you’ll need a GitHub Personal Access Token (PAT):

* Go to: https://github.com/settings/tokens
* Click **Generate new token (classic)**
* Select **repo** scope
* Copy the token

Then in R:

```r
Sys.setenv(GITHUB_PAT = "your_token_here")
```

### 3. Install the package

```r
remotes::install_github("uri-ncipher/NetClusterSim")
```

### 4. Load the package

```r
library(NetClusterSim)
```

---

