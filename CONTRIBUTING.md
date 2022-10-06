# Contributing

## Releasing

1. Commit your changes:

    ```bash
    git commit -a
    ```

2. Choose a version number e.g.:

    ```bash
    version=release-1.1.0
    ```

3. Update the [changelog](./CHANGELOG.md):

    ```bash
    git commit -a -m "Update changelog for $version"
    ```

4. Tag the release:

    ```bash
    git tag -a $version -m "release $version"
    ```

5. Push the changes to the GitHub repository

    ```bash
    git push -u origin main --follow-tags
    ```