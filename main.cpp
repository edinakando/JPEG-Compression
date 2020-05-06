#include "stdafx.h"
#include "common.h"
#include <vector>

using namespace std;

#define BLOCK_SIZE 8

int quantizationVals[] = { 16, 11, 10, 16, 24, 40, 51, 61,
						  12, 12, 14, 19, 26, 58, 60, 55,
						  14, 13, 16, 24, 40, 57, 69, 56,
						  14, 17, 22, 29, 51, 87, 80, 62,
						  18, 22, 37, 56, 68, 109, 103, 77,
						  24, 35, 55, 64, 81, 104, 113, 92,
						  49, 64, 78, 87, 103, 121, 120, 101,
						  72, 92, 95, 98, 112, 100, 103, 99 };
Mat_<int> quantization(BLOCK_SIZE, BLOCK_SIZE, quantizationVals);

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

	for (int i = 0; i < img.rows; i += BLOCK_SIZE)
		for (int j = 0; j < img.cols; j += BLOCK_SIZE) {
			cv::Rect rectangle = cv::Rect(i, j, BLOCK_SIZE, BLOCK_SIZE);
			pixels.push_back(cv::Mat_<Vec3b>(img, rectangle));
		}

	return pixels;
}

double alpha(int i) {
	if (i == 0)
		return sqrt(1.0 / 8);
	return sqrt(2.0 / 8);
}

Mat_<double> applyFDCT(Mat_<int> img) {
	Mat_<double> C(img.rows, img.cols);
	img -= 128; //center around 0

	for(int i = 0; i < img.rows; i++)
		for (int j = 0; j < img.cols; j++) {
			double sum = 0;
			for (double x = 0; x < BLOCK_SIZE; x++)
				for (double y = 0; y < BLOCK_SIZE; y++)
					sum += img(x, y) * cos(PI * (2 * x + 1) * i / (2 * BLOCK_SIZE))
									 * cos(PI * (2 * y + 1) * j / (2 * BLOCK_SIZE));

			C(i, j) = alpha(i) * alpha(j) * sum;
		}

	return C;
}

Mat_<uchar> applyIDCT(Mat_<int> img) {
	Mat_<int> V(img.rows, img.cols);

	for(int i = 0; i < img.rows; i++)
		for (int j = 0; j < img.cols; j++) {
			double sum = 0;
			for (double x = 0; x < BLOCK_SIZE; x++)
				for (double y = 0; y < BLOCK_SIZE; y++)
					sum += img(x, y) * cos(PI * (2 * y + 1) * i / (2 * BLOCK_SIZE))
									 * cos(PI * (2 * x + 1) * j / (2 * BLOCK_SIZE));

			V(i, j) = round(alpha(i) * alpha(j) * sum);
		}

	V += 128;
	Mat_<uchar> result = V.clone();
	return result;
}

Mat_<int> applyQuantization(Mat_<double> src) {
	Mat_<int> dst(src.rows, src.cols);

	for(int i = 0; i < src.rows; i++)
		for (int j = 0; j < src.cols; j++) 
			dst(i, j) = round(src(i, j) / quantization(i, j));

	return dst;
}


Mat_<int> applyDequantization(Mat_<int> src) {
	Mat_<int> dst(src.rows, src.cols);

	for (int i = 0; i < src.rows; i++)
		for (int j = 0; j < src.cols; j++)
			dst(i, j) = src(i, j) * quantization(i, j);

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

vector<pair<int, int>> runLengthEncode(vector<int> vals) {
	vector<pair<int, int>> res;
	int i = 0;
	while (i < vals.size()) {
		int zeroesCount = 0;
		while (i < vals.size() - 1 && vals[i] == 0) {
			zeroesCount++;
			i++;
		}

		if (i == vals.size() - 1 && vals[i++] == 0)   //the rest is 0
			res.push_back(pair<int, int>(0, 0));
		else res.push_back(pair<int, int>(zeroesCount, vals[i++]));
	}

	return res;
}

void huffmanCode() {

}

vector<vector<pair<int, int>>> compressImage(Mat_<Vec3b> src) {
	Mat_<uchar> Y(src.rows, src.cols), Cr(src.rows, src.cols), Cb(src.rows, src.cols);
	Mat_<Vec3b> dst(src.rows, src.cols);

	cvtColor(src, dst, COLOR_BGR2YCrCb);
	vector<Mat> channels;
	split(dst, channels);
	Y = channels[0];
	Cr = channels[1];
	Cb = channels[2];

	vector<Mat_<int>> yBlocks = getDivisionInBlocks(Y);

	vector<vector<pair<int, int>>> RLEs;
	for (Mat_<uchar> block : yBlocks) {
		Mat_<double> dct = applyFDCT(block);
		Mat_<int> quantized = applyQuantization(dct);
		vector<int> zigZag = zigZagTraversal(quantized);
		RLEs.push_back(runLengthEncode(zigZag));
	}

	//todo: huffman?

	return RLEs;
}

Mat_<uchar> decompress(Mat_<int> quantization) {
	//todo: get quantized val from huffman encoding
	Mat_<int> deq = applyDequantization(quantization);
	Mat_<uchar> decompressedImage = applyIDCT(deq);
	return decompressedImage;
}


int main()
{
	/*Mat_<Vec3b> src = openImage();
	compressImage(src);*/

	/*int vals[] = { 64, 56, 56, 57, 70, 84, 84, 59,
				  66, 64, 35, 36, 87, 45, 21, 58,
				  66, 66, 66, 59, 35, 87, 26, 104,
				  35, 75, 76, 45, 81, 37, 34, 35,
				  45, 96, 125, 107, 31, 15, 107, 90,
				  88, 89, 88, 78, 64, 57, 85, 81,
				  62, 59, 68, 113, 144, 104, 66, 73,
				  107, 121, 89, 21, 35, 64, 65, 65 };*/

	
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
	cout << "Original:\n" << Y << endl << endl;
	cout <<  "Normalized around 0:\n " << Y << endl << endl;

	Mat_<double> rez = applyFDCT(Y);
	cout << "FDCT:\n" << rez << endl << endl;

	Mat_<int> q = applyQuantization(rez);
	cout << "Quantized:\n" << q << endl << endl;

	//vector<int> rr = zigZagTraversal(q);
	//cout << "ZigZag: \n";
	//for (int i = 0; i < rr.size(); i++)
	//	cout << rr[i] << " ";
	//cout << endl << endl;

	//cout << "RLE:\n";
	//vector<pair<int, int>> res1 = runLengthEncode(rr);
	//for (int i = 0; i < res1.size(); i++)
	//	cout << res1[i].first << ", " << res1[i].second << endl;
	//cout << endl << endl;

	Mat_<int> deq = applyDequantization(q);
	cout << "Dequantized: \n"  << deq << endl << endl;

	Mat_<double> ddd = applyIDCT(deq);
	cout << "IDCT:\n" <<  ddd << endl << endl;

	//imshow("original", src);
	//imshow("Y", Y);
	//imshow("Cr", Cr);
	//imshow("Cb", Cb);

	//waitKey();
}