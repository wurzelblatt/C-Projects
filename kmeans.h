struct Kmeans{
	//Solve for a k-means clustering model from a set of data points and initial guesses of the means.
	//Output is a set of means and an assignment of each data point to one component.
	Int nn, mm, kk, nchg;
	MatDoub data, means;
	VecInt assign, count;
	Kmeans(MatDoub &data, MatDoub &mmeans) : nn(ddata.nrows()), mm(ddata.ncols()),
	kk(means.nrows()), data(ddata), means(mmeans), assign(nn), count(kk) {
		//Constructor. Arguments are the data points (as rows in a matrix), and initial guesses for
		//the means (also as rows in a matrix).
		estep();														//Perform one initial E-step and M-step.																
		mstep();														//User is responsible for calling additional
																		//steps until convergence is obtained.
	}
	Int estep() {
		//Perform one E-step.
		Int k,m,n,kmin;
		Doub dmin,d;
		nchg = 0.;
		for 
