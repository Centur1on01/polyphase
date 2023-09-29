#define _USE_MATH_DEFINES

#include <iostream>
#include <complex>
#include <vector>
#include <deque>
#include <ctime>
#include <cmath>
#include <random>
#include <string>
#include <fstream>
#include <cstdint>

using namespace std;



class signalSource {
	// наследуемые классы: источник бита и маппер
	// источник бита берет из файла или m-последовательности
	// мапппер позволяет выбрать сигнальное созвездие

	// получение битов из файла -> выбор созвездия -> перевод битов в complex
	//string fileName;
	double Ps;
	complex<double> result;

public:

	signalSource(const double& Ps = 2) : Ps(Ps) {
		//result = { sqrt(Ps / 2.), 0.0 };
		result = { Ps, 0.0 };
		srand(time(0));
	}

	complex<double> next() {
		return (rand() % 2 == 0 ? result : -result);
	}

};



class upsampling {
	// вставка нулей 2/4/8
	int n;
	int upsampling_n;
	complex<double> signalup;

public:

	upsampling(const int& upsampling_n = 1) : upsampling_n(upsampling_n) {
		n = 0;
	}

	complex<double> next(const complex<double>& signal) {
		if (n % upsampling_n == 0)
			signalup = signal;
		else
			signalup = 0;
		++n;
		return signalup;
	}

	int getUpsampling_n() {
		return upsampling_n;
	}

};



class FIR {
	// прием коэффициентов из файла
	// функция next для фильтрации
	double order;
	double coeff;
	complex<double> filteredSignal, temp;
	vector<double> coeffTable;
	deque<complex<double>> buff;
	string fileName;
	ifstream file;

public:

	FIR(const string& fileName = "filterTable.txt") : fileName(fileName), order(0), filteredSignal(0) {
		ifstream file(fileName, ios::in);

		if (file.is_open())
			while (!file.eof()) {
				file >> coeff;
				/*if (coeff == NULL)
					break;*/
				coeffTable.push_back(coeff);
				++order;
			}
		//--order;
		file.close();
		buff.resize(order, 0);
	}

	// цикл со сверткой с отдельной переменной которая отсчитывает свертку для коэфов
	// возвращает отфильтрованный отсчет
	complex<double> next(const complex<double>& signal) {
		buff.push_back(signal);
		buff.pop_front();

		filteredSignal = { 0.0, 0.0 }; // 
		for (int i = 0; i < buff.size(); ++i) {
			//filteredSignal += coeffTable[i] * buff[buff.size() - i - 1]; // coeffTable[i] * buff[buff.size() - i - 1];
			filteredSignal += coeffTable[i] * buff[i]; // coeffTable[i] * buff[buff.size() - i - 1];
		}
		//temp = filteredSignal;

		return filteredSignal;
	}

	double getOrder() {
		return order;
	}

	void clear() {
		buff.resize(order, 0);
	}

};



class polyphaseInterpolation {
	double order;
	double coeff;
	double interpolationCoeff;
	double multiplyNumber;
	double it, j;
	complex<double> filteredSignal, temp;
	vector<double> coeffTable;
	vector<complex<double>> sum;
	deque<complex<double>> buff, signalbuff, decimbuff;
	string fileName;
	ifstream file;

public:

	polyphaseInterpolation(const string& fileName = "filterTable.txt", double interpolationCoeff = 2) : fileName(fileName), interpolationCoeff(interpolationCoeff), order(0), filteredSignal(0) {
		ifstream file(fileName, ios::in);

		if (file.is_open())
			while (!file.eof()) {
				file >> coeff;
				coeffTable.push_back(coeff);
				++order;
			}
		file.close();
		buff.resize(order, 0);
		//multiplyNumber = order / interpolationCoeff;
		//sum.resize(order);
		//signalbuff.resize(multiplyNumber, 0);
		//signalbuff.resize(interpolationCoeff);
		it = 0;
		j = 0;
		//sum.resize(order / interpolationCoeff, 0);
	}

	void filter(deque<complex<double>> signal) {
		signalbuff.resize(interpolationCoeff, 0);
		for (int i = 0; i < interpolationCoeff; ++i) {
			filteredSignal = 0;
			for (int j = 0; j < order / interpolationCoeff; ++j) {
				filteredSignal += signal[j] * coeffTable[order - 1 - (interpolationCoeff * (j)) - i];
				
			}
			signalbuff.push_back(filteredSignal);
			signalbuff.pop_front();
		}
	}

	complex<double> getFilteredSignal() {
		complex<double> removed = signalbuff.back();
		signalbuff.pop_back();
		return removed;
	}

	double getOrder() {
		return order;
	}

	//complex<double> filt(complex<double> signal) {
	//	if (j == decimationCoeff - 1) {
	//		j = 0;
	//	}
	//	for (int i = 0; i < decimationCoeff; ++i) {
	//		decimbuff[i] += signal * coeffTable[order - 1 - j - decimationCoeff * i];
	//	}


	//}

	//complex<double> interpolationFilter(complex<double> signal) {
	//	/*for (int i = 0; i < (order - 1) / interpolationCoeff; ++i) {
	//		for (int j = 0; i < interpolationCoeff; ++j) {

	//		}
	//	}*/

	//	if (it < multiplyNumber) {
	//		signalbuff.push_back(signal);
	//		signalbuff.pop_front();
	//		++it;
	//		return 0;
	//	}
	//	else {

	//	}


	//	for (int i = 0; i < multiplyNumber; ++i) {
	//		for (int j = 0; j < interpolationCoeff; ++j) {
	//			sum[i * 2 + j] += 1; 
	//		}
	//	}

	//}

	//complex<double> getInterpolatedSignal() {
	//	complex<double> temp = sum[0];
	//	buff.push_back(0);
	//	buff.pop_front();
	//	return temp;
	//}




};



class polyphaseDecimation {
	double order;
	double coeff;
	double decimationCoeff;
	double multiplyNumber;
	double it, j;
	complex<double> filteredSignal, temp;
	vector<double> coeffTable;
	vector<complex<double>> sum;
	deque<complex<double>> buff, signalbuff, decimbuff;
	string fileName;
	ifstream file;

public:

	polyphaseDecimation(const string& fileName = "filterTable.txt", double decimationCoeff = 2) : fileName(fileName), decimationCoeff(decimationCoeff), order(0), filteredSignal(0) {
		ifstream file(fileName, ios::in);

		if (file.is_open())
			while (!file.eof()) {
				file >> coeff;
				coeffTable.push_back(coeff);
				++order;
			}
		file.close();
		buff.resize(order, 0);
		//multiplyNumber = order / interpolationCoeff;
		//sum.resize(order);
		//signalbuff.resize(multiplyNumber, 0);
		signalbuff.resize(decimationCoeff);
		it = 0;
		j = 0;
		//sum.resize(order / interpolationCoeff, 0);
	}

	complex<double> filter(deque<complex<double>> signal) {
		filteredSignal = 0;
		for (int i = 0; i < decimationCoeff; ++i) {
			for (int j = 0; j < order / decimationCoeff; ++j) {
				filteredSignal += signal[i] * coeffTable[order - 1 - (decimationCoeff * (j)) - i];
			}
		}
		


		return filteredSignal;
	}

};



// измерение мощности сигнала
class PWR {
	double L, result, temp;
	deque<complex<double>> buff;

public:

	PWR(const double& L = 64) : L(L), result(0) {
		buff.resize(L, 0);
	}

	// возвращает сумму квадратов амплитуд отсчетов, деленную на L
	double next(const complex<double>& signal) {
		//if (buff.size() != L)
		//	buff.push_front(signal); // 
		//else {
		//	buff.pop_back(); //
		//	buff.push_front(signal); //
		//}

		buff.push_back(signal);
		result -= pow(buff[0].real(), 2) + pow(buff[0].imag(), 2);
		buff.pop_front();

		for (int i = 0; i < buff.size(); ++i) {
			result += pow(buff[L - i - 1].real(), 2) + pow(buff[L - i - 1].imag(), 2);
		}
		temp = result / L;
		result = 0;

		return temp;
	}

	double getPWR() {
		return temp;
	}


};



class awgn {
	// генерация шума next
	// pwr measure класс и сам шум
	default_random_engine generator;
	normal_distribution<double> distr;
	double Pn, upsampling_n;

public:

	awgn(const double& Pn = 2, const double& upsampling_n = 1) : Pn(Pn), upsampling_n(upsampling_n) {
		this->Pn *= upsampling_n;
		distr = normal_distribution<double>(0.0, 1.0);
	}

	double next(bool measured = 0, double snr = 0, double PWR = 0) {
		if (measured == 0)
			return std::sqrt(Pn / upsampling_n / 2.) * distr(generator);
		else {
			Pn = pow(10., (10. * log10(PWR) - snr) / 10.) * upsampling_n; // upsampling_n не нужен, потому что PWR делится на upsampling_n в формирующем фильтре
			return std::sqrt(Pn / 2.) * distr(generator);
		}
	}

};



class decision {
	// >= 0 = 1, < 0 = -1
	complex<double> result;
	double Ps;

public:

	decision(const double& Ps = 2) : Ps(Ps) {}

	// возвращает демодулированный отсчёт
	complex<double> next(const complex<double>& signal) {
		if (signal.real() >= 0.0)
			return { Ps, 0.0 };
		else
			return { -Ps, 0.0 };
	}

};



class biterror {
	// сравнение decision выхода и изначального signal через функцию next
	double N, correct;

public:

	biterror(const double& N = 1) : N(N), correct(0) {}

	// возвращает значение BER
	double next(complex<double> original, complex<double> final) {
		if (original == final)
			++correct;
		return (N - correct) / N;
	}

	void clear() {
		correct = 0;
	}

};



int main()
{

	int N = 1000;
	//int upsampling_n = 2;
	int L = 64;

	double Ps = 1.;
	double PsdB = 10. * log10(Ps);
	//double POWER = Ps / upsampling_n;
	double Pn, PndB, BERRatio = 0, BERRatio_t = 0;



	complex<double> signal, signalProcessed, signalFinal;

	signalSource source(Ps);
	//upsampling upsampler(upsampling_n);
	//FIR shapingFilter("filterTable_2.txt");
	//FIR matchedFilter("filterTable_2.txt");
	PWR measurer(L);
	decision decisionmaker(Ps);
	biterror biterr(N);
	//biterror biterr_t(N);

	//double delay = shapingFilter.getOrder();
	//double delayPerFilter = (delay - 1) / upsampling_n;
	deque<complex<double>> buff_orig, buff_fin;

	//ofstream test("testAfter1.txt");
	ofstream out("out.txt");
	//ofstream snrcheck("snr.txt");
	//ofstream pwrout("pwr.txt");
	ofstream berout("ber.txt");

	ofstream PCMM("PCM.pcm", ios::binary);
	//ofstream matlab("matlab.txt");
	//ofstream demodcheck("demod.txt");
	//ofstream demodcheck2("demod2.txt");

	double it = 0;

	deque<complex<double>> signalAccumulated;
	int interpolationCoeff = 4;
	int decimationCoeff = interpolationCoeff;
	polyphaseInterpolation pintf("intTable4_160.txt", interpolationCoeff); // intTable2_160.txt, intTable4_160.txt, intTable8_160.txt
	polyphaseDecimation pdecf("intTable4_160.txt", decimationCoeff);
	int interpolationFilterOrder = pintf.getOrder();
	deque<complex<double>> signalInterpolated, signalPrevious;
	deque<complex<double>> checkk, signalDelayed;
	complex<double> signalDecimated, signalDecided;

	//bool ready = 0;
	//bool firstElementDone = 0;

	//checkk.resize(interpolationFilterOrder / 4 - 28, 0); // order / 4 | interpolationCoeff + 2
	//signalDelayed.resize(decimationCoeff, 0);



	double POWER = Ps / (double)interpolationCoeff;

	for (int snr = 0; snr < 1; ++snr) {
		PndB = (PsdB - snr) - 11;
		Pn = pow(10., PndB / 10.);
		awgn noise(Pn, interpolationCoeff);

		signalAccumulated.resize(0);
		signalInterpolated.resize(interpolationCoeff, 0);
		signalPrevious.resize(interpolationCoeff, 0);
		/*
			intcoeff = 2 : checkk.resize(interpolationFilterOrder / 4 + 2, 0);
			intcoeff = 4 : checkk.resize(interpolationFilterOrder / 4 - 18, 0);
			intcoeff = 8 : checkk.resize(interpolationFilterOrder / 4 - 28, 0);
		*/
		checkk.resize(interpolationFilterOrder / 4 - 18, 0);
		signalDelayed.resize(decimationCoeff, 0);

		double accumulated = 0;
		bool ready = 0;
		bool firstElementDone = 0;

		for (int i = 0; i < N; ++i) {
			signal = source.next();
			checkk.push_back(signal);
			checkk.pop_front();
			++accumulated;

			if (!ready)
				signalAccumulated.push_back(signal);

			if (accumulated == interpolationFilterOrder / interpolationCoeff) {
				signalAccumulated.pop_front();
				signalAccumulated.push_back(signal);
				//pintf.filter(signalAccumulated);
				ready = 1;
			}

			if (ready) {

				pintf.filter(signalAccumulated);
				signalAccumulated.pop_front();
				signalAccumulated.push_back(signal);

				for (int j = 0; j < interpolationCoeff; ++j) {
					complex<double> tmp = pintf.getFilteredSignal();

					// шум
					//tmp = { tmp.real() + noise.next(0), tmp.imag() + noise.next(0) };
					//tmp = { tmp.real() + noise.next(1, snr, POWER), tmp.imag() + noise.next(1, snr, POWER) };

					/*float reald = (float)tmp.real();
					float imagd = (float)tmp.imag();
					PCMM.write((const char*)&reald, sizeof(float));
					PCMM.write((const char*)&imagd, sizeof(float));*/

					signalInterpolated.push_back(tmp);
					signalInterpolated.pop_front();
				}

				if (!firstElementDone) {
					signalDelayed.pop_front();
					signalDelayed.push_back(signalInterpolated[0]);
					firstElementDone = 1;
					//signalPrevious = signalInterpolated;
				}
				else {
					//signalDelayed.pop_front();
					//signalDelayed.push_back(signalInterpolated[signalInterpolated.size()]);
					for (int k = 0; k < decimationCoeff - 1; ++k) { // < decimationCoeff - 1
						signalDelayed.pop_front();
						signalDelayed.push_back(signalPrevious[k + 1]); // signalInterpolated | signalPrevious || decimationCoeff - 1 - k
					}
					signalDelayed.pop_front();
					signalDelayed.push_back(signalInterpolated[0]);
				}

				signalPrevious = signalInterpolated;

				signalDecimated = pdecf.filter(signalDelayed); // signalInterpolated

				signalDecided = decisionmaker.next(signalDecimated);

				out << i << "\t" << signal << "\t" << signalDecided << "\t" << checkk[0] << endl;

				//BERRatio = biterr.next(signal, signalDecided);
				BERRatio = biterr.next(checkk[0], signalDecided);
			}
		}

		berout << snr << "\t" << BERRatio << "\t" << 0.5 * erfc(sqrt(pow(10., snr / 10.))) << endl;


		biterr.clear();
	}

	int aaa = 1; // stop


	out.close();
	//snrcheck.close();
	//pwrout.close();
	berout.close();
	PCMM.close();
	//demodcheck.close();
	//demodcheck2.close();

	return 0;
}
