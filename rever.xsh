$PROJECT = 'diffpy.structure'
$ACTIVITIES = [
              'tag',  # Creates a tag for the new version number
              'push_tag',  # Pushes the tag up to the $TAG_REMOTE
              'pypi',  # Sends the package to pypi
              'ghrelease'  # Creates a Github release entry for the new tag
               ]
$PUSH_TAG_REMOTE = 'git@github.com:diffpy/diffpy.structure.git'  # Repo to push tags to
$GITHUB_ORG = 'diffpy'  # Github org for Github releases and conda-forge
$GITHUB_REPO = 'diffpy.structure'  # Github repo for Github releases  and conda-forge
