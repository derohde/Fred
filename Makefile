pre:
	sudo apt install -y libboost-all-dev
	sudo apt-get install -y python3-setuptools
	sudo apt-get install -y python3-numpy
	sudo apt-get install -y python-setuptools
	sudo apt-get install -y python-numpy
	sudo apt-get install -y cmake

python3:
	cd py && python3 ./setup.py install --user
	
python2:
	cd py && python2 ./setup.py install --user

clean:
	pip3 uninstall Fred -y
	pip2 uninstall Fred -y
	rm -r py/dist py/build/ py/Fred.egg-info/
