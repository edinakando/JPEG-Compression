#include "stdafx.h"
#include "common.h"
#include <vector>
#include <fstream>

using namespace std;

#define BLOCK_SIZE 8

int imgRows, imgCols;

int quantizationVals[] = { 16, 11, 10, 16, 24, 40, 51, 61,
						  12, 12, 14, 19, 26, 58, 60, 55,
						  14, 13, 16, 24, 40, 57, 69, 56,
						  14, 17, 22, 29, 51, 87, 80, 62,
						  18, 22, 37, 56, 68, 109, 103, 77,
						  24, 35, 55, 64, 81, 104, 113, 92,
						  49, 64, 78, 87, 103, 121, 120, 101,
						  72, 92, 95, 98, 112, 100, 103, 99 };
Mat_<int> quantization(BLOCK_SIZE, BLOCK_SIZE, quantizationVals);

bool isInside(Mat img, int i, int j) {
	return i >= 0 && i < img.rows&& j >= 0 && j < img.cols;
}

Mat openImage() {
	char fname[MAX_PATH];
	openFileDlg(fname);
	Mat src = imread(fname, CV_LOAD_IMAGE_COLOR);
	waitKey();
	return src;
}

Mat_<Vec3b> padImage(Mat_<Vec3b> src) {
	int rows = src.rows;
	int cols = src.cols;
	while (rows % BLOCK_SIZE != 0) rows++;
	while (cols % BLOCK_SIZE != 0) cols++;
	Mat_<Vec3b> dst(rows, cols);

	for (int i = 0; i < src.rows; i++)
		for (int j = 0; j < src.cols; j++)
			dst(i, j) = src(i, j);

	for (int i = src.rows; i < rows; i++)
		for (int j = src.cols; j < cols; j++)
			dst(i, j) = 0;

	return dst;
}

vector<Mat_<int>> getDivisionInBlocks(Mat_<int> img) {
	vector<Mat_<int>> pixels;

	for (int i = 0; i < img.rows; i += BLOCK_SIZE)
		for (int j = 0; j < img.cols; j += BLOCK_SIZE) {
			Mat_<int> block(BLOCK_SIZE, BLOCK_SIZE);

			for (int u = 0; u < BLOCK_SIZE; u++) {
				int currentI = i + u;
				for (int v = 0; v < BLOCK_SIZE; v++) {
					int currentJ = j + v;

					if (isInside(img, currentI, currentJ))
						block(u, v) = img(currentI, currentJ);
				}
			}
			pixels.push_back(block);
		}

	return pixels;
}

double alpha(int i) {
	if (i == 0)
		return sqrt(1.0 / 8);
	return sqrt(2.0 / 8);
}

Mat_<int> applyFDCT(Mat_<int> img) {
	Mat_<int> DCT(img.rows, img.cols);
	img -= 128; //center around 0

	for (int i = 0; i < img.rows; i++)
		for (int j = 0; j < img.cols; j++) {
			double sum = 0;
			for (double x = 0; x < BLOCK_SIZE; x++)
				for (double y = 0; y < BLOCK_SIZE; y++)
					sum += img(x, y) * cos(PI * (2 * x + 1) * i / (2 * BLOCK_SIZE))
					* cos(PI * (2 * y + 1) * j / (2 * BLOCK_SIZE));

			DCT(i, j) = alpha(i) * alpha(j) * sum;
		}

	return DCT;
}

Mat_<uchar> applyIDCT(Mat_<int> DCT) {
	Mat_<int> IDCT(DCT.rows, DCT.cols, { 0 });

	for (int i = 0; i < DCT.rows; i++)
		for (int j = 0; j < DCT.cols; j++) {
			for (int x = 0; x < BLOCK_SIZE; x++)
				for (int y = 0; y < BLOCK_SIZE; y++)
					IDCT(i, j) += alpha(x) * alpha(y) * DCT(x, y) 
								* cos(PI * (2 * i + 1) * x / (2 * BLOCK_SIZE))
								* cos(PI * (2 * j + 1) * y / (2 * BLOCK_SIZE));
		}

	IDCT += 128;
	return IDCT;
}

Mat_<int> applyQuantization(Mat_<int> src, Mat_<int> quantizationTable) {
	Mat_<int> dst(src.rows, src.cols);

	for (int i = 0; i < src.rows; i++)
		for (int j = 0; j < src.cols; j++)
			dst(i, j) = round((float)src(i, j) / quantizationTable(i, j));

	return dst;
}

Mat_<int> applyDequantization(Mat_<int> src, Mat_<int> quantizationTable) {
	Mat_<int> dst(src.rows, src.cols);

	for (int i = 0; i < src.rows; i++)
		for (int j = 0; j < src.cols; j++)
			dst(i, j) = round((float)src(i, j) * quantizationTable(i, j));

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

Mat_<int> reconstructZigZagTraversal(vector<int> traversal) {
	Mat_<int> original(BLOCK_SIZE, BLOCK_SIZE);

	int currentIndex = 0;
	for (int i = 0; i < 2 * BLOCK_SIZE - 1; i++) {
		if (i % 2 == 0) {
			int currentRow = i < BLOCK_SIZE ? i : BLOCK_SIZE - 1;
			int currentColumn = i < BLOCK_SIZE ? 0 : i - BLOCK_SIZE + 1;

			while (currentRow >= 0 && currentColumn < BLOCK_SIZE)
				original(currentRow--, currentColumn++) = traversal[currentIndex++];
		}
		else {
			int currentRow = i < BLOCK_SIZE ? 0 : i - BLOCK_SIZE + 1;
			int currentColumn = i < BLOCK_SIZE ? i : BLOCK_SIZE - 1;

			while (currentColumn >= 0 && currentRow < BLOCK_SIZE && i < traversal.size())
				original(currentRow++, currentColumn--) = traversal[currentIndex++];
		}
	}

	return original;
}

vector<int> runLengthEncode(vector<int> vals) {
	vector<int> rle;
	int i = 0;
	while (i < vals.size()) {
		int zeroesCount = 0;
		while (i < vals.size() && vals[i] == 0) {
			zeroesCount++;
			i++;
		}

		if (zeroesCount != 0) {
			rle.push_back(0);
			rle.push_back(zeroesCount);
		}
		else {
			rle.push_back(vals[i]);
			i++;
		}
	}

	return rle;
}

vector<int> decodeRunLength(vector<int> rle) {
	vector<int> decodedRLE;
	int i = 0;

	while (i < rle.size()) {
		if (rle[i] == 0) {
			decodedRLE.insert(decodedRLE.end(), rle[i + 1], rle[i]);
			i++;
		}
		else decodedRLE.push_back(rle[i]);

		i++;
	}

	return decodedRLE;
}

vector<vector<int>> compressImageComponent(Mat_<int> component, Mat_<int> quantizationTable) {
	vector<Mat_<int>> blocks = getDivisionInBlocks(component);

	vector<vector<int>> RLEs;
	for (Mat_<uchar> block : blocks) {
		Mat_<int> dct = applyFDCT(block);
		Mat_<int> quantized = applyQuantization(dct, quantizationTable);
		vector<int> zigZag = zigZagTraversal(quantized);
		RLEs.push_back(runLengthEncode(zigZag));
	}

	return RLEs;
}

Mat_<uchar> decompressBlock(vector<int> rle, Mat_<int> quantizationTable) {
	vector<int> zigZagTraversal = decodeRunLength(rle);
	Mat_<int> quantization = reconstructZigZagTraversal(zigZagTraversal);
	Mat_<int> deq = applyDequantization(quantization, quantizationTable);
	return applyIDCT(deq);
}

Mat_<uchar> reconstructImageComponent(vector<vector<int>> rles, Mat_<int> quantizationTable) {
	Mat_<uchar> img(imgRows, imgCols);

	int currentBlockIndex = 0;
	for (int i = 0; i < img.rows; i += BLOCK_SIZE)
		for (int j = 0; j < img.cols; j += BLOCK_SIZE) {
			Mat_<uchar> currentBlock = decompressBlock(rles[currentBlockIndex++], quantizationTable);

			for (int u = 0; u < BLOCK_SIZE; u++)
				for (int v = 0; v < BLOCK_SIZE; v++) {
					int currentI = i + u;
					int currentJ = j + v;

					if (isInside(img, currentI, currentJ))
						img(currentI, currentJ) = currentBlock(u, v);
				}
		}

	return img;
}

int main()
{
	Mat_<Vec3b> src = openImage();
	imgRows = src.rows;
	imgCols = src.cols;

	int originalSize = src.rows * src.cols;
	int compressedSize = 0;

	Mat_<Vec3b> padded = padImage(src);

	Mat_<uchar> Y(padded.rows, padded.cols), Cr(padded.rows, padded.cols), Cb(padded.rows, padded.cols);
	Mat_<Vec3b> dst(padded.rows, padded.cols);

	cvtColor(padded, dst, COLOR_BGR2YCrCb);
	Mat channels[3];
	split(dst, channels);
	Y = channels[0];
	Cr = channels[1];
	Cb = channels[2];

	vector<vector<int>> yRLEs = compressImageComponent(Y, 2 * quantization);
	vector<vector<int>> crRLEs = compressImageComponent(Cr, 4 * quantization);
	vector<vector<int>> cbRLEs = compressImageComponent(Cb, 4 * quantization);

#pragma region FILE_WRITE
	ofstream file;
	file.open("img.txt");
	file << padded.rows << " " << padded.cols << " ";

	for (int i = 0; i < yRLEs.size(); i++) {
		for (int j = 0; j < yRLEs[i].size(); j++) {
			file << yRLEs[i][j] << " ";
			compressedSize++;
		}
		file << endl;
	}
	for (int i = 0; i < crRLEs.size(); i++) {
		for (int j = 0; j < crRLEs[i].size(); j++) {
			file << crRLEs[i][j] << " ";
			compressedSize++;
		}
		file << endl;
	}
	for (int i = 0; i < cbRLEs.size(); i++) {
		for (int j = 0; j < cbRLEs[i].size(); j++) {
			file << cbRLEs[i][j] << " ";
			compressedSize++;
		}
		file << endl;
	}
	file.close();
#pragma endregion

	cout << "Original Size: " << originalSize << endl;
	cout << "Compressed Size: " << compressedSize << endl;
	cout << "Compression Ratio: " << (float) originalSize / compressedSize << endl;

#pragma region FILE_READ
	ifstream input;
	input.open("img.txt");

	int rows, cols;
	input >> rows >> cols;

	int blockNo = (rows / BLOCK_SIZE) * (cols / BLOCK_SIZE);
	
	vector<vector<int>> yRead, crRead, cbRead;
	int x;
	string  line;

	for (int i = 0; i < blockNo; i++) {
		vector<int> current;
		getline(input, line);
		stringstream linestream(line);
		while (linestream >> x)
			current.push_back(x);
		yRead.push_back(current);
	}
	for (int i = 0; i < blockNo; i++) {
		vector<int> current;
		getline(input, line);
		stringstream linestream(line);
		while (linestream >> x)
			current.push_back(x);
		crRead.push_back(current);
	}
	for (int i = 0; i < blockNo; i++) {
		vector<int> current;
		getline(input, line);
		stringstream linestream(line);
		while (linestream >> x)
			current.push_back(x);
		cbRead.push_back(current);
	}

#pragma endregion

	Mat_<uchar> yReconstructed = reconstructImageComponent(yRead, 2 * quantization);
	Mat_<uchar> crReconstructed = reconstructImageComponent(crRead, 4 * quantization);
	Mat_<uchar> cbReconstructed = reconstructImageComponent(cbRead, 4 * quantization);

	channels[0] = yReconstructed;
	channels[1] = crReconstructed;
	channels[2] = cbReconstructed;

	Mat mergedRes, finalRes;
	merge(channels, 3, mergedRes);
	cvtColor(mergedRes, finalRes, COLOR_YCrCb2BGR);

	imshow("original", src);
	imshow("reconstructed", finalRes);

	waitKey();
}