#include "stdafx.h"
#include "common.h"
#include <vector>

using namespace std;

#define PIXEL_SIZE 8

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

vector<Mat_<Vec3b>> getDivisionInPixels(Mat_<uchar> img) {
	vector<Mat_<Vec3b>> pixels;
	for (int i = 0; i < img.rows; i += PIXEL_SIZE)
		for (int j = 0; j < img.cols; j += PIXEL_SIZE) {
			cv::Rect rectangle = cv::Rect(i, j, PIXEL_SIZE, PIXEL_SIZE);
			pixels.push_back(cv::Mat_<Vec3b>(img, rectangle));
		}

	return pixels;
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

	vector<Mat_<Vec3b>> yPixels = getDivisionInPixels(Y);
	vector<Mat_<Vec3b>> crPixels = getDivisionInPixels(Cr);
	vector<Mat_<Vec3b>> cbPixels = getDivisionInPixels(Cb);

	//TODO - dct
	imshow("original", src);
	imshow("Y", Y);
	imshow("Cr", Cr);
	imshow("Cb", Cb);

	waitKey();
}