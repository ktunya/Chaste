class StepSizeException : public std::exception {
	
public:

	double 	displacement;
	double 	suggestedNewStep;
	const std::string message;
	bool isTerminal;

	StepSizeException(double errDisplacement, double errNewStep, const std::string errMessage, bool errIsTerminal):
		displacement(errDisplacement),
		suggestedNewStep(errNewStep),
		message(errMessage),
		isTerminal(errIsTerminal),
		std::exception()
	{
	}

	virtual const char* what() const throw() {
		return message.c_str();
	}

	~StepSizeException() throw() {} 
};