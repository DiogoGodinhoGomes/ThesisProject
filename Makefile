calc: func.py data.py
	@pep8 --show-source --show-pep8 func.py
	@./func.py

	@pep8 --show-source --show-pep8 data.py
	@./data.py

show:
	@ds9 -multiframe -zscale Images/*.fits &
