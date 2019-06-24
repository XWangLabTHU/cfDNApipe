## Python package: pysam, numpy, os, sys, collections
echo "===========================================================================>>"
echo "                                                                             "
echo "Detecting python package 'numpy', 'pysam'......"
echo "                                                                             "
if python -c 'import pkgutil; exit(not pkgutil.find_loader("numpy"))'; then
	echo 'Package numpy exist!'
	flag_numpy=1
else
	echo 'Package numpy not exist!'
	flag_numpy=0
fi

if python -c 'import pkgutil; exit(not pkgutil.find_loader("pysam"))'; then
	echo 'Package pysam exist!'
	flag_pysam=1
else
	echo 'Package pysam not exist!'
	flag_pysam=0
fi

if [[ $flag_numpy -eq 1 && $flag_pysam -eq 1 ]]; then
	echo "Python packages exist!"
else
	echo "Some packages not exist."
	read -p "Do you want to install them (y or n): " flag_input
	if [ "$flag_input" = "n" ]; then
		echo "Some python packages are not available, exit now!"
		exit 1
	elif [ "$flag_input" = "y" ]; then
		if [ $flag_numpy -eq 0 ]; then
			echo "                                                                             "
			echo "Now, installing numpy......"
			conda install -c anaconda numpy 2>&1
			flag_install=$?
			if [ $flag_install -eq 0 ]; then
				echo "Numpy installation complete!"
			else
				echo "Numpy installation failed, exit now."
				exit 1
			fi
		fi
		if [ $flag_pysam -eq 0 ]; then
			echo "                                                                             "
			echo "Now, installing pysam......"
			conda install -c bioconda pysam 2>&1
			flag_install=$?
			if [ $flag_install -eq 0 ]; then
				echo "Pysam installation complete!"
			else
				echo "Pysam installation failed, exit now."
				exit 1
			fi
		fi
	else
		echo "Input must be y or n!"
		exit 1
	fi
fi
echo "                                                                             "
echo "<<==========================================================================="