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

vector<Mat_<uchar>> getDivisionInBlocks(Mat_<uchar> img) {
	vector<Mat_<uchar>> pixels;
	Mat_<uchar> normalizedImg = img.clone();
	normalizedImg -= 128;
	for (int i = 0; i < normalizedImg.rows; i += BLOCK_SIZE)
		for (int j = 0; j < normalizedImg.cols; j += BLOCK_SIZE) {
			cv::Rect rectangle = cv::Rect(i, j, BLOCK_SIZE, BLOCK_SIZE);
			pixels.push_back(cv::Mat_<Vec3b>(normalizedImg, rectangle));
		}

	return pixels;
}

double alpha(int n) {
	if (n == 0)
		return sqrt(1 / n);
	return sqrt(2 / n);
}

Mat_<double> convertToFrequencyUsingFDCT(Mat_<uchar> img) {
	Mat_<double> C(img.rows, img.cols);

	for(int u = 0; u < img.rows; u++)
		for (int v = 0; v < img.cols; v++) {
			double sum = 0;
			for (double x = 0; x < BLOCK_SIZE; x++)
				for (double y = 0; y < BLOCK_SIZE; y++)
					sum += img(x, y) * cos(PI * (2 * x + 1) * u / (2 * BLOCK_SIZE))
									 * cos(PI * (2 * x + 1) * v / (2 * BLOCK_SIZE));

			C(u, v) = alpha(u) * alpha(v) * sum;
		}

	return C;
}

int main()
{
	Mat_<Vec3b> src = openImage();
	Mat_<Vec3b> dst(src.rows, src.cols);
	Mat_<uchar> Y(src.rows, src.cols), Cr(src.rows, src.cols), Cb(src.rows, src.cols);

	cvtColor(src, dst, COLOR_BGR2YCrCb);
	vector<Mat> channels;
	split(dst, channels);
	Y = channels[0];
	Cr = channels[1];
	Cb = channels[2];

	vector<Mat_<uchar>> yBlocks = getDivisionInBlocks(Y);
	vector<Mat_<uchar>> crBlocks = getDivisionInBlocks(Cr);
	vector<Mat_<uchar>> cbBlocks = getDivisionInBlocks(Cb);

	vector<Mat_<double>> dct;
	for (Mat_<uchar> block : yBlocks)
		dct.push_back(convertToFrequencyUsingFDCT(block));

	//TODO - apply to cr and cb?

	//TODO - dct
	imshow("original", src);
	imshow("Y", Y);
	imshow("Cr", Cr);
	imshow("Cb", Cb);

	waitKey();
}