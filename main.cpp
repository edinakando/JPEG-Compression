#include "stdafx.h"
#include "common.h"
#include <vector>

using namespace std;

#define BLOCK_SIZE 8

Mat openImage() {
	char fname[MAX_PATH];
	openFileDlg(fname);
	Mat src = imread(fname, CV_LOAD_IMAGE_COLOR);
	waitKey();
	return src;
}

void encodePixels(vector<Mat_<Vec3b>> pixels) {
	for (Mat_<Vec3b> pixel : pixels) {

	}
}

vector<Mat_<int>> getDivisionInBlocks(Mat_<uchar> img) {
	vector<Mat_<int>> pixels;
	Mat_<int> normalizedImg = img.clone();
	normalizedImg -= 128; 	//center around 0
	for (int i = 0; i < normalizedImg.rows; i += BLOCK_SIZE)
		for (int j = 0; j < normalizedImg.cols; j += BLOCK_SIZE) {
			cv::Rect rectangle = cv::Rect(i, j, BLOCK_SIZE, BLOCK_SIZE);
			pixels.push_back(cv::Mat_<Vec3b>(normalizedImg, rectangle));
		}

	return pixels;
}

double alpha(int i) {
	if (i == 0)
		return sqrt(1.0 / 8);
	return sqrt(2.0 / 8);
}

Mat_<double> convertToFrequencyUsingFDCT(Mat_<int> img) {
	Mat_<double> C(img.rows, img.cols);

	for(int u = 0; u < img.rows; u++)
		for (int v = 0; v < img.cols; v++) {
			double sum = 0;
			for (double x = 0; x < BLOCK_SIZE; x++)
				for (double y = 0; y < BLOCK_SIZE; y++)
					sum += img(x, y) * cos(PI * (2 * x + 1) * u / (2 * BLOCK_SIZE))
									 * cos(PI * (2 * y + 1) * v / (2 * BLOCK_SIZE));

			C(u, v) = alpha(u) * alpha(v) * sum;
			double a = alpha(u);
			double b = C(u, v);
			int c = 5;
		}

	return C;
}

Mat_<int> applyQuantization(Mat_<double> src, Mat_<int> quantization) {
	Mat_<int> dst(src.rows, src.cols);

	for(int i = 0; i < src.rows; i++)
		for (int j = 0; j < src.cols; j++) 
			dst(i, j) = round(src(i, j) / quantization(i, j));

	return dst;
}

vector<int> zigZagTraversal(Mat_<int> src) {
	vector<int> res;

	for (int i = 0; i < src.rows + src.cols - 1; i++) {
		if (i % 2 == 0) {
			int currentRow = i < src.rows ? i : src.rows - 1;
			int currentColumn = i < src.cols ? 0 : i - src.cols + 1;

			while (currentRow >= 0 && currentColumn < src.cols)
				res.push_back(src(currentRow--, currentColumn++));
		}
		else {
			int currentRow = i < src.rows ? 0 : i - src.rows + 1;
			int currentColumn = i < src.cols ? i : src.cols - 1;

			while (currentColumn >= 0 && currentRow < src.rows) 
				res.push_back(src(currentRow++, currentColumn--));
		}
	}

	return res;
}

int main()
{

	/*int vals[] = { 64, 56, 56, 57, 70, 84, 84, 59,
				  66, 64, 35, 36, 87, 45, 21, 58,
				  66, 66, 66, 59, 35, 87, 26, 104,
				  35, 75, 76, 45, 81, 37, 34, 35,
				  45, 96, 125, 107, 31, 15, 107, 90,
				  88, 89, 88, 78, 64, 57, 85, 81,
				  62, 59, 68, 113, 144, 104, 66, 73,
				  107, 121, 89, 21, 35, 64, 65, 65 };*/

	//Mat_<Vec3b> src = openImage();
	//Mat_<Vec3b> dst(src.rows, src.cols);
	//Mat_<uchar> Y(src.rows, src.cols), Cr(src.rows, src.cols), Cb(src.rows, src.cols);

	//cvtColor(src, dst, COLOR_BGR2YCrCb);
	//vector<Mat> channels;
	//split(dst, channels);
	//Y = channels[0];
	//Cr = channels[1];
	//Cb = channels[2];

	//vector<Mat_<int>> yBlocks = getDivisionInBlocks(Y);
	//vector<Mat_<int>> crBlocks = getDivisionInBlocks(Cr);
	//vector<Mat_<int>> cbBlocks = getDivisionInBlocks(Cb);

	//vector<Mat_<double>> dct;
	//for (Mat_<uchar> block : yBlocks)
	//	dct.push_back(convertToFrequencyUsingFDCT(block));

	////TODO - apply to cr and cb?


	//testing
	int vals[] = { 52, 55, 61, 66, 70, 61, 64, 73,
				  64, 59, 55, 90, 109, 85, 69, 72,
				  62, 59, 68, 113, 144, 104, 66, 73,
				  63, 58, 71, 122, 154, 106, 70, 69,
				  67, 61, 68, 104, 126, 88, 68, 70,
				  79, 65, 60, 70, 77, 68, 58, 75,
				  85, 71, 64, 59, 55, 61, 65, 83,
				  87, 79, 69, 68, 65, 76, 78, 94 };
	Mat_<int> Y(8, 8, vals);
	Y -= 128;
	cout << Y << endl;

	Mat_<double> rez = convertToFrequencyUsingFDCT(Y);
	cout << rez;
	int quantizationVals[] = {16, 11, 10, 16, 24, 40, 51, 61,
							  12, 12, 14, 19, 26, 58, 60, 55,
							  14, 13, 16, 24, 40, 57, 69, 56,
							  14, 17, 22, 29, 51, 87, 80, 62,
							  18, 22, 37, 56, 68, 109, 103, 77,
							  24, 35, 55, 64, 81, 104, 113, 92,
							  49, 64, 78, 87, 103, 121, 120, 101,
							  72, 92, 95, 98, 112, 100, 103, 99};
	Mat_<int> quantization(BLOCK_SIZE, BLOCK_SIZE, quantizationVals);
	Mat_<int> q = applyQuantization(rez, quantization);

	cout << endl;
	cout << q << endl;

	vector<int> rr =  zigZagTraversal(q);
	for (int i = 0; i < rr.size(); i++)
		cout << rr[i] << " ";

	////TODO - dct
	//imshow("original", src);
	//imshow("Y", Y);
	//imshow("Cr", Cr);
	//imshow("Cb", Cb);

	//waitKey();
}