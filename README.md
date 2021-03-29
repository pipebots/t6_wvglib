```
 __________.__             ___.     |__|  __
 \______   \__|_____   ____\_ |__   _||__/  |_  ______ (C) George Jackson-Mills 2020
  |     ___/  \____ \_/ __ \| __ \ /  _ \   __\/  ___/
  |    |   |  |  |_> >  ___/| \_\ (  O_O )  |  \___ \
  |____|   |__|   __/ \___  >___  /\____/|__| /____  >
              |__|        \/    \/                 \/
```

# wvglib - Suite of functions for the analysis of waveguides

## Overview

This library contains the Python implementation of the equations used to analyse waveguides. Currently the focus is on lossy waveguides, i.e. ones embedded in a dielectric. There is also a rudimentary module on circular waveguides with a conductive boundary. Plans for the future include add rectangular metal waveguides.

## Requirements

As my current `Python` development environment is a bit polluted, I'll list the required packages here. Apologies.

- `Python>=3.6`
- `numpy`
- `scipy`
- `rflib` -> in-house developed Python package, see repo

## Tests

To be added. Bad practice, I know.

## Installation

Use `pip install -e .` in the folder to which you clone or download this. This will install `rflib` as an "editable" package in your current environment, meaning you should just do a `git pull` in the future to get any updates.

## Contributing

Contributions are more than welcome and are in fact actively sought! Please contact Viktor either at [eenvdo@leeds.ac.uk](mailto:eenvdo@leeds.ac.uk) or on the `Pipebots` Slack.