# New M4 Software Installation Guide

### M4

Firstly, clone the official M4 software from `GitHub` (updated fequently with respect to GitLab):

```bash
https://github.com/ArcetriAdaptiveOptics/M4.git
```

Inside the clone repo (assuming it's in `~/git/M4`), switch to the most recent branch, which is:

```bash
git checkout pf-SW-cleanup
```

At this point, install the package by:

```bash
pip install --use-pep517 .
```

### Dependencies (OptiCalib)

For the `OptiCalib` dependency, it is more straightforward. Just:

```bash
pip install git+/https://github.com/pietroferraiuolo/labott.git
```

## Software startup

### opticalib

To start the working environment for M4, firstly you need to set up the opticalib package. With the package, a script has been installed, called `calpy` (check the info by running `calpy -h` in your terminal). Choose a path where you want to store your configurations and data. Let's assume it will be `~/experiments/m4`.
Then:

```bash
calpy -c experiments/m4
```

This command will create the folder structure, as well as the configuration file needed for the software (in `~/experiments/m4/SysConfig/configuration.yaml`).

Finally, to initiate a working session with everythin set up, run

```bash
calpy -f experiments/m4
```

This command will start an ipython session, with all the relevant paths and configurations loaded into the opticalib package. Inside the shell, the following modules of the `opticalib` package are available:

- `opticalib` (alias `opt`)
- `dmutils`
- `zernike` (alias `zern`)
- `osutils` (alias `osu`)

_NOTE_: The two above steps can be done compressed into a single command:

```bash
calpy -f experiments/m4 --create
```

This command will create the folders and configuration file, and start the ipython session.

### m4

After setting up the opticalib package, with it's configuration file

```python
import m4
```

Importing m4 will modify the `configuration.,yaml` file, to add a flag for simulated (or not) devices, in order to make the simulator work (It will be soon removed...)

... And everything is set.

The available functions/classes in the M4 package (as is) are:

- `create_ott` - functions that returns the `ott`, the `dm` and the `interferometer` instances
- `OttAligner` - A wrapper of the more general `opticalib.alignment.Alignment` class, specific for the OTT
- `OttIffAcquisition` - A wrapper for iff acquisition and analysis specific for the DP. _Work in Progress..._

## Docs

For a comprehensive documentation about the configuration file, please check [this link](https://github.com/pietroferraiuolo/labott/blob/main/opticalib/core/_configurations/DOCS.md).