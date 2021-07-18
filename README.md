# Atmospheric Turubulence Simulator

Simulates atmospheric turbulence for astronomical observations.

## Download

This is a module of the [milk package](https://github.com/milk-org/milk).

To download source code :
~~~
git clone https://github.com/milk-org/milk.git
cd milk
./fetch_atmosphere_dev.sh
~~~

## Compile

From within milk source code root directory :

~~~
mkdir _build
cd _build
cmake ..
make
sudo make install
~~~

## Generating atmospheric wavefronts

From user work directory.

Create configuration file :
~~~
milk-atmturb-mkconf WFsim.conf
~~~

Edit file as needed, then run computation for desired wavelength value (unit micron) :

~~~
milk-atmturb-runturb 1.65
~~~

The process will run until killed (CTRL-C) and write wavefronts to disk.

