# chandra123

# DESCRIPTION:
A set of scripts for downloading and analyzing Chandra/ACIS image files for detecting or placing upper limits on the presence of emission at a set of input coordinates.

# REQUIREMENTS:
In addition to the requirements in requirements.txt, this script requires the latest version (4.11) of ciao (to use ciao/make_instmap_weights and merge_obs, specifically).

ciao can be installed from: http://cxc.harvard.edu/ciao/threads/ciao_install_tool/

After installing ciao, you can set up the remainder of the chandra123 environment with the following commands:

```
alias ciao='source $CIAO_DIR/ciao-4.11/bin/ciao.bash'
ciao
pip3 freeze > $ASCDS_INSTALL/constraints.txt
pip3 install -c $ASCDS_INSTALL/constraints.txt 'astropy<3.1' requests xmltodict
```
