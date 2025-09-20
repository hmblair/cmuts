## Installation

### 1. Obtain `cmuts`

Clone the `cmuts` repository with all submodules:

```bash
git clone --recurse-submodules https://github.com/hmblair/cmuts
cd cmuts
```


### 2. Build `cmuts`

The configuration script contains all required setup and build steps.

```
./configure
```

If building the multithreaded version, pass the `--mpi` flag to the script.

```
./configure --mpi
```

A debug build, which reduces optimization and enables various sanitizers, is also possible.

```
./configure --debug
```

!!! warning
    Building with `--mpi` or `--debug` after an existing build has already been performed will require deleting the `build` directory first.

### 3. Modify PATH

Afterwards, add the `cmuts` binary directory to your PATH:

=== "Temporary (current session)"
    ```bash
    export PATH="$(pwd)/bin:$PATH"
    ```

=== "Permanent (bash)"
    ```bash
    echo 'export PATH="'$(pwd)'/bin:$PATH"' >> ~/.bashrc
    source ~/.bashrc
    ```

=== "Permanent (zsh)"
    ```bash
    echo 'export PATH="'$(pwd)'/bin:$PATH"' >> ~/.zshrc
    source ~/.zshrc
    ```

## Updating

To update to a newer version, from within the `cmuts` directory, run

```bash
git pull origin master
git submodule update --recursive
rm -rf build
./configure # OR ./configure --mpi
```
