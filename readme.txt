Integrated RNA-seq Analysis Application (IRNAA)
================
Roseric AZONDEKON

March 11, 2019
--------------------


INSTALLATION
----------------
1. Extract IRNAA to your local drive.

2. To install IRNAA on Ubuntu
--------------------------------
 From the command line interface, change your directory (using the `cd` command) to the extracted local `irnaa` folder on your local drive, and run the following script:

sudo bash linux_install.sh && source ~/.bashrc

You should now launch the IRNAA shiny application by executing the `irnaa` command from your shell terminal:

irnaa


2. To install IRNAA on Windows
--------------------------------

before installing IRNAA, make sure that R is installed and can be run from the command line.

To make R accessible from the command line interface, you may use the following instructions:

- Search for the Rscript.exe executable and copy its location path (e.g C:\Program Files\R\R-3.5.1\bin)
- Open the start menu and type in "View advanced system settings", click on "Environment variables"
- Under "System variables", select Path and click on edit
    - Windows 10 and related Click "New", and paste in the directory path to Rscript.exe.
    - For earlier versions of Windows systems, paste in the directory path at the beginning of the path line followed by a semicolon (no space before or after the semicolon)
- Open the Windows Command Prompt and run Rscript --version which should output the version of R installed on your computer.

Now, from the extracted local irnaa folder, right-click on the win_install.bat file and run it as administrator. When the installation is succesful, a shortcut to IRNAA is created on your desktop.

To launch IRNAA, double-click on the IRNAA shortcut on your desktop.

To install it from the DOS command line interface, open your command prompt (as administrator), change your directory (using the cd command) to the extracted local irnaa folder and execute the following code:

win_install

You may also run IRNAA by double-clicking on the irnaa.bat file located the win directory inside your locally extracted irnaa folder.
