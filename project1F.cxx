#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <math.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>

using std::cerr;
using std::endl;

#define NORMALS

double ceil_441(double f)
{
    return ceil(f-0.00001);
}

double floor_441(double f)
{
    return floor(f+0.00001);
}

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
   std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInputData(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}

struct LightingParameters
{
	LightingParameters(void)
	{
		lightDir[0] = -0.6;
		lightDir[1] = 0;
		lightDir[2] = -0.8;
		Ka = 0.3;
		Kd = 0.7;
		Ks = 2.3;
		alpha = 2.5;
	};


	double lightDir[3]; // The direction of the light source
	double Ka;           // The coefficient for ambient lighting.
	double Kd;           // The coefficient for diffuse lighting.
	double Ks;           // The coefficient for specular lighting.
	double alpha;        // The exponent term for specular lighting.
};

LightingParameters lp;

class Triangle
{
  public:
      double         X[3];
      double         Y[3];
	  double		 Z[3];
	  double		 colors[3][3];
	  double		 normals[3][3];

// would some methods for transforming the triangle in place be helpful?
};

class Screen
{
public:
	unsigned char   *buffer;
	int width, height;

	// would some methods for accessing and setting pixels be helpful?
};

std::vector<Triangle>
GetTriangles(void)
{
	vtkPolyDataReader *rdr = vtkPolyDataReader::New();
	rdr->SetFileName("proj1e_geometry.vtk");
	//cerr << "Reading" << endl;
	rdr->Update();
	//cerr << "Done reading" << endl;
	if (rdr->GetOutput()->GetNumberOfCells() == 0)
	{
		cerr << "Unable to open file!!" << endl;
		exit(EXIT_FAILURE);
	}
	vtkPolyData *pd = rdr->GetOutput();

	int numTris = pd->GetNumberOfCells();
	vtkPoints *pts = pd->GetPoints();
	vtkCellArray *cells = pd->GetPolys();
	vtkDoubleArray *var = (vtkDoubleArray *)pd->GetPointData()->GetArray("hardyglobal");
	double *color_ptr = var->GetPointer(0);
	//vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
	//float *color_ptr = var->GetPointer(0);
	vtkFloatArray *n = (vtkFloatArray *)pd->GetPointData()->GetNormals();
	float *normals = n->GetPointer(0);
	std::vector<Triangle> tris(numTris);
	vtkIdType npts;
	vtkIdType *ptIds;
	int idx;
	for (idx = 0, cells->InitTraversal(); cells->GetNextCell(npts, ptIds); idx++)
	{
		if (npts != 3)
		{
			cerr << "Non-triangles!! ???" << endl;
			exit(EXIT_FAILURE);
		}
		double *pt = NULL;
		pt = pts->GetPoint(ptIds[0]);
		tris[idx].X[0] = pt[0];
		tris[idx].Y[0] = pt[1];
		tris[idx].Z[0] = pt[2];
#ifdef NORMALS
		tris[idx].normals[0][0] = normals[3 * ptIds[0] + 0];
		tris[idx].normals[0][1] = normals[3 * ptIds[0] + 1];
		tris[idx].normals[0][2] = normals[3 * ptIds[0] + 2];
#endif
		pt = pts->GetPoint(ptIds[1]);
		tris[idx].X[1] = pt[0];
		tris[idx].Y[1] = pt[1];
		tris[idx].Z[1] = pt[2];
#ifdef NORMALS
		tris[idx].normals[1][0] = normals[3 * ptIds[1] + 0];
		tris[idx].normals[1][1] = normals[3 * ptIds[1] + 1];
		tris[idx].normals[1][2] = normals[3 * ptIds[1] + 2];
#endif
		pt = pts->GetPoint(ptIds[2]);
		tris[idx].X[2] = pt[0];
		tris[idx].Y[2] = pt[1];
		tris[idx].Z[2] = pt[2];
#ifdef NORMALS
		tris[idx].normals[2][0] = normals[3 * ptIds[2] + 0];
		tris[idx].normals[2][1] = normals[3 * ptIds[2] + 1];
		tris[idx].normals[2][2] = normals[3 * ptIds[2] + 2];
#endif

		// 1->2 interpolate between light blue, dark blue
		// 2->2.5 interpolate between dark blue, cyan
		// 2.5->3 interpolate between cyan, green
		// 3->3.5 interpolate between green, yellow
		// 3.5->4 interpolate between yellow, orange
		// 4->5 interpolate between orange, brick
		// 5->6 interpolate between brick, salmon
		double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
		double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
		unsigned char RGB[8][3] = { { 71, 71, 219 },
									{ 0, 0, 91 },
									{ 0, 255, 255 },
									{ 0, 128, 0 },
									{ 255, 255, 0 },
									{ 255, 96, 0 },
									{ 107, 0, 0 },
									{ 224, 76, 76 }
		};
		for (int j = 0; j < 3; j++)
		{
			float val = color_ptr[ptIds[j]];
			int r;
			for (r = 0; r < 7; r++)
			{
				if (mins[r] <= val && val < maxs[r])
					break;
			}
			if (r == 7)
			{
				cerr << "Could not interpolate color for " << val << endl;
				exit(EXIT_FAILURE);
			}
			double proportion = (val - mins[r]) / (maxs[r] - mins[r]);
			tris[idx].colors[j][0] = (RGB[r][0] + proportion * (RGB[r + 1][0] - RGB[r][0])) / 255.0;
			tris[idx].colors[j][1] = (RGB[r][1] + proportion * (RGB[r + 1][1] - RGB[r][1])) / 255.0;
			tris[idx].colors[j][2] = (RGB[r][2] + proportion * (RGB[r + 1][2] - RGB[r][2])) / 255.0;
		}
	}

	return tris;
}

double dotProduct(double v1[3], double v2[3])
{
	double product = 0;
	for (int i = 0; i < 3; ++i)
	{
		product += v1[i] * v2[i];
	}

	return product;
}

double *crossProduct(double v1[3], double v2[3])
{
	double *retv = (double *)malloc(sizeof(double) * 3);
	retv[0] = v1[1] * v2[2] - v1[2] * v2[1];
	retv[1] = v2[0] * v1[2] - v2[2] * v1[0];
	retv[2] = v1[0] * v2[1] - v1[1] * v2[0];

	return retv;
}

double *normalize(double *v1)
{
	double norm = sqrt((v1[0] * v1[0]) + (v1[1] * v1[1]) + (v1[2] * v1[2]));
	double *retv = (double *)malloc(sizeof(double) * 3);
	retv[0] = v1[0] / norm;
	retv[1] = v1[1] / norm;
	retv[2] = v1[2] / norm;

	return retv;
}

class Matrix
{
public:
	double          A[4][4];  // A[i][j] means row i, column j

	void            TransformPoint(const double *ptIn, double *ptOut);
	static Matrix   ComposeMatrices(const Matrix &, const Matrix &);
	void            Print(ostream &o);
};

void
Matrix::Print(ostream &o)
{
	for (int i = 0; i < 4; i++)
	{
		char str[256];
		sprintf(str, "(%.7f %.7f %.7f %.7f)\n", A[i][0], A[i][1], A[i][2], A[i][3]);
		o << str;
	}
}

Matrix
Matrix::ComposeMatrices(const Matrix &M1, const Matrix &M2)
{
	Matrix rv;
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
		{
			rv.A[i][j] = 0;
			for (int k = 0; k < 4; k++)
				rv.A[i][j] += M1.A[i][k] * M2.A[k][j];
		}

	return rv;
}

void
Matrix::TransformPoint(const double *ptIn, double *ptOut)
{
	ptOut[0] = ptIn[0] * A[0][0]
		+ ptIn[1] * A[1][0]
		+ ptIn[2] * A[2][0]
		+ ptIn[3] * A[3][0];
	ptOut[1] = ptIn[0] * A[0][1]
		+ ptIn[1] * A[1][1]
		+ ptIn[2] * A[2][1]
		+ ptIn[3] * A[3][1];
	ptOut[2] = ptIn[0] * A[0][2]
		+ ptIn[1] * A[1][2]
		+ ptIn[2] * A[2][2]
		+ ptIn[3] * A[3][2];
	ptOut[3] = ptIn[0] * A[0][3]
		+ ptIn[1] * A[1][3]
		+ ptIn[2] * A[2][3]
		+ ptIn[3] * A[3][3];
}

double *subtract(double *v1, double*v2)
{
	double *retv = (double *)malloc(sizeof(double) * 3);
	retv[0] = v1[0] - v2[0];
	retv[1] = v1[1] - v2[1];
	retv[2] = v1[2] - v2[2];
	return retv;
}

class Camera
{
public:
	double          near, far;
	double          angle;
	double          position[3];
	double          focus[3];
	double          up[3];

	Matrix          ViewTransform(void);
	Matrix          CameraTransform(void);
	Matrix          DeviceTransform(double, double);
};

Matrix
Camera::ViewTransform()
{
	Matrix m;
	m.A[0][0] = 1/tan(angle/2);
	m.A[0][1] = 0;
	m.A[0][2] = 0;
	m.A[0][3] = 0;

	m.A[1][0] = 0;
	m.A[1][1] = 1 / tan(angle / 2);
	m.A[1][2] = 0;
	m.A[1][3] = 0;

	m.A[2][0] = 0;
	m.A[2][1] = 0;
	m.A[2][2] = (far + near) / (far - near);
	m.A[2][3] = -1;

	m.A[3][0] = 0;
	m.A[3][1] = 0;
	m.A[3][2] = (2*far*near)/(far-near);
	m.A[3][3] = 0;

	return m;
}

Matrix
Camera::CameraTransform()
{
	double *O = position;
	double *u = normalize(crossProduct(up, subtract(O,focus)));
	double *v = normalize(crossProduct(subtract(O, focus), u));
	double *w = normalize(subtract(O,focus));
	double t[3];
	t[0] = 0 - O[0];
	t[1] = 0 - O[1];
	t[2] = 0 - O[2];

	//printf("U = %f, %f, %f\n", u[0], u[1], u[2]);
	//printf("V = %f, %f, %f\n", v[0], v[1], v[2]);
	//printf("W = %f, %f, %f\n", w[0], w[1], w[2]);

	Matrix m;
	m.A[0][0] = u[0];
	m.A[0][1] = v[0];
	m.A[0][2] = w[0];
	m.A[0][3] = 0;

	m.A[1][0] = u[1];
	m.A[1][1] = v[1];
	m.A[1][2] = w[1];
	m.A[1][3] = 0;

	m.A[2][0] = u[2];
	m.A[2][1] = v[2];
	m.A[2][2] = w[2];
	m.A[2][3] = 0;

	m.A[3][0] = dotProduct(u, t);
	m.A[3][1] = dotProduct(v, t);
	m.A[3][2] = dotProduct(w, t);
	m.A[3][3] = 1;

	free(u);
	free(v);
	free(w);

	return m;
}

Matrix
Camera::DeviceTransform(double n, double m)
{
	Matrix j;
	j.A[0][0] = n / 2;
	j.A[0][1] = 0;
	j.A[0][2] = 0;
	j.A[0][3] = 0;

	j.A[1][0] = 0;
	j.A[1][1] = m / 2;
	j.A[1][2] = 0;
	j.A[1][3] = 0;

	j.A[2][0] = 0;
	j.A[2][1] = 0;
	j.A[2][2] = 1;
	j.A[2][3] = 0;

	j.A[3][0] = n / 2;
	j.A[3][1] = m / 2;
	j.A[3][2] = 0;
	j.A[3][3] = 1;

	return j;

}

double SineParameterize(int curFrame, int nFrames, int ramp)
{
	int nNonRamp = nFrames - 2 * ramp;
	double height = 1. / (nNonRamp + 4 * ramp / M_PI);
	if (curFrame < ramp)
	{
		double factor = 2 * height*ramp / M_PI;
		double eval = cos(M_PI / 2 * ((double)curFrame) / ramp);
		return (1. - eval)*factor;
	}
	else if (curFrame > nFrames - ramp)
	{
		int amount_left = nFrames - curFrame;
		double factor = 2 * height*ramp / M_PI;
		double eval = cos(M_PI / 2 * ((double)amount_left / ramp));
		return 1. - (1 - eval)*factor;
	}
	double amount_in_quad = ((double)curFrame - ramp);
	double quad_part = amount_in_quad * height;
	double curve_part = height * (2 * ramp) / M_PI;
	return quad_part + curve_part;
}

Camera
GetCamera(int frame, int nframes)
{
	double t = SineParameterize(frame, nframes, nframes / 10);
	Camera c;
	c.near = 5;
	c.far = 200;
	c.angle = M_PI / 6;
	c.position[0] = 40 * sin(2 * M_PI*t);
	c.position[1] = 40 * cos(2 * M_PI*t);
	c.position[2] = 40;
	c.focus[0] = 0;
	c.focus[1] = 0;
	c.focus[2] = 0;
	c.up[0] = 0;
	c.up[1] = 1;
	c.up[2] = 0;
	return c;
}


double **OrderPoints(double *p1, double *p2, double *p3, double colors[3][3], double shading[3])
{
	double **points = (double**)malloc(sizeof(double*) * 7);
	double newshade[3];
	if (p1[1] == p2[1])
	{
		points[2] = p3;
		points[3] = colors[2];
		if (p1[0] < p2[0])
		{
			points[0] = p1;
			points[1] = colors[0];
			points[4] = p2;
			points[5] = colors[1];
			newshade[0] = shading[0];
			newshade[1] = shading[2];
			newshade[2] = shading[1];
			points[6] = newshade;
		}
		else
		{
			points[0] = p2;
			points[1] = colors[1];
			points[4] = p1;
			points[5] = colors[0];
			newshade[0] = shading[1];
			newshade[1] = shading[2];
			newshade[2] = shading[0];
			points[6] = newshade;
		}
	}
	else if (p2[1] == p3[1])
	{
		points[2] = p1;
		points[3] = colors[0];
		if (p2[0] < p3[0])
		{
			points[0] = p2;
			points[1] = colors[1];
			points[4] = p3;
			points[5] = colors[2];
			newshade[0] = shading[1];
			newshade[1] = shading[0];
			newshade[2] = shading[2];
			points[6] = newshade;
		}
		else
		{
			points[0] = p3;
			points[1] = colors[2];
			points[4] = p2;
			points[5] = colors[1];
			newshade[0] = shading[2];
			newshade[1] = shading[0];
			newshade[2] = shading[1];
			points[6] = newshade;
		}
	}
	else if (p1[0] < p3[0])
	{
		points[0] = p1;
		points[1] = colors[0];
		points[2] = p2;
		points[3] = colors[1];
		points[4] = p3;
		points[5] = colors[2];
		newshade[0] = shading[0];
		newshade[1] = shading[1];
		newshade[2] = shading[2];
		points[6] = newshade;
	}
	else
	{
		points[0] = p3;
		points[1] = colors[2];
		points[2] = p2;
		points[3] = colors[1];
		points[4] = p1;
		points[5] = colors[0];
		newshade[0] = shading[2];
		newshade[1] = shading[1];
		newshade[2] = shading[0];
		points[6] = newshade;
	}
	return points;
}

double slope(double *pt1, double *pt2)
{
	if (pt1[0] == pt2[0])
		return 0.0;
	return ((pt2[1] - pt1[1]) / (pt2[0] - pt1[0]));
}

void colorPixel(unsigned char *buffer, double *colors, int width, int x, int y, double z, double *zbuffer, double shading)
{
	double r, g, b;
	//check if pixels are trying to overlap onto left side
	if (x == 1000 && y >= 0)
		return;
	//check depth buffer to see if we should display these
	if (z >= zbuffer[(width*y) + x])
	{
		zbuffer[(width*y) + x] = z;

		r = (ceil_441(colors[0] * 255 * shading));
		g = (ceil_441(colors[1] * 255 * shading));
		b = (ceil_441(colors[2] * 255 * shading));

		if (r > 255)
			r = 255;
		if (g > 255)
			g = 255;
		if (b > 255)
			b = 255;
		buffer[(3 * width*y) + 3 * x] = r;
		buffer[(3 * width*y) + 3 * x + 1] = g;
		buffer[(3 * width*y) + 3 * x + 2] = b;
	}
}

double *interpolateColors(double *p1, double *p1color, double *p2, double *p2color, double end)
{
	double *newcolors = (double*)malloc(sizeof(double) * 3);
	double t;

	t = (end - p1[1]) / (p2[1] - p1[1]);

	newcolors[0] = (p1color[0] + (t*(p2color[0] - p1color[0])));
	newcolors[1] = (p1color[1] + (t*(p2color[1] - p1color[1])));
	newcolors[2] = (p1color[2] + (t*(p2color[2] - p1color[2])));
	return newcolors;
}

//Return the type of trianle based on the vertexes
char *TriType(double *p1, double *p2, double *p3)
{
	if (p1[1] == p2[1])
	{
		if (p3[1] > p1[1])
			return "goingup";
		else return "goingdown";
	}
	else if (p2[1] == p3[1])
	{
		if (p1[1] > p2[1])
			return "goingup";
		else return "goingdown";
	}
	else if (p1[1] == p3[1])
	{
		if (p2[1] > p1[1])
			return "goingup";
		else return "goingdown";
	}
	else return "arbitrary";
}

//Calculate pixels for any going down triangles
void RasterizeGoingDown(double *p1, double *p2, double *p3, unsigned char *buffer, double *zbuffer, double width, double color[3][3], double shading[3])
{
//printf("Inside Going Down Rasterize\n");
	double **points;
	double *topl, *topr, *bottom;
	double *toplcolor, *toprcolor, *bottomcolor;
	double lslope, rslope, lend, rend, lb, rb, lshade, rshade, shade;
	double *lcolors, *rcolors;
	double *midcolors = (double*)malloc(sizeof(double) * 3);
	double lt, rt, t, lz, rz, z;
	double toplshade, toprshade, bottomshade;

	//return points in order followed by color
	points = OrderPoints(p1, p2, p3, color, shading);
	topl = points[0];
	toplcolor = points[1];
	toplshade = points[6][0];

	bottom = points[2];
	bottomcolor = points[3];
	bottomshade = points[6][1];

	topr = points[4];
	toprcolor = points[5];
	toprshade = points[6][2];

	double rowmin = ceil_441(bottom[1]);
	double rowmax = floor_441(topl[1]);
	//make sure triangles are in bounds
	if (rowmax >= 1000)
		rowmax = 999;
	if (rowmin < 0)
		rowmin = 0;

	lslope = slope(topl, bottom);
	rslope = slope(bottom, topr);

//printf("ROWMIN: %f, ROWMAX: %f\n", rowmin, rowmax);

	for (double j = rowmin; j <= rowmax; j++)
	{
		if (lslope == 0.0)
			lend = topl[0];
		else
		{
			lb = (topl[1] - (lslope*topl[0]));
			lend = (j - lb) / lslope;
		}
		//check possibiility of out of bounds wherever possible
		if (lend < 0.0)
			lend = 0.0;
		if (lend > 1000)
			lend = 1000;
		if (rslope == 0.0)
			rend = bottom[0];
		else
		{
			rb = (topr[1] - (rslope*topr[0]));
			rend = (j - rb) / rslope;
		}
		if (rend > 1000)
			rend = 1000;
		if (rend < 0)
			rend = 0;

		//Interpolate colors using y values from points
		lcolors = interpolateColors(topl, toplcolor, bottom, bottomcolor, j);
		rcolors = interpolateColors(bottom, bottomcolor, topr, toprcolor, j);

		//Interpolate the z value for the leftend
		lt = (j - topl[1]) / (bottom[1] - topl[1]);

		lz = topl[2] + (lt*(bottom[2] - topl[2]));
		lshade = toplshade + (lt * (bottomshade - toplshade));

		//Interpolate the z value for the rightend
		rt = (j - bottom[1]) / (topr[1] - bottom[1]);

		rz = bottom[2] + (rt*(topr[2] - bottom[2]));
		rshade = bottomshade + (rt * (toprshade - bottomshade));
		
//printf("lz = %f + (%f * (%f - %f))\n", topl[2], lt, bottom[2], topl[2]);
//printf("rz = %f + (%f * (%f - %f))\n", bottom[2], rt, topr[2], bottom[2]);

		//Go through all points along the line
		for (double k = ceil_441(lend); k <= floor_441(rend); k++)
		{
			//calculate colors for the individual pixel
			t = (k - lend) / (rend - lend);
			midcolors[0] = lcolors[0] + (t*(rcolors[0] - lcolors[0]));
			midcolors[1] = lcolors[1] + (t*(rcolors[1] - lcolors[1]));
			midcolors[2] = lcolors[2] + (t*(rcolors[2] - lcolors[2]));

			//calculate depth to pass to colorPixel
			z = lz + (t * (rz - lz));
			shade = (lshade + (t *(rshade - lshade)));

			//color pixel
			colorPixel(buffer, midcolors, width, k, j, z, zbuffer, shade);
		}
		free(lcolors);
		free(rcolors);
	}
	free(midcolors);
	free(points);
}
//Calculate pixels for any going up triangles
void RasterizeGoingUp(double *p1, double *p2, double *p3, unsigned char *buffer, double *zbuffer, double width, double color[3][3], double shading[3])
{
//printf("Inside Going UP Rasterize\n");
	double **points;
	double *botl, *botr, *top;
	double lslope, rslope, lend, rend, lb, rb, t;
	double *botlcolor, *topcolor, *botrcolor;
	double *lcolors, *rcolors;
	double lt, rt, lz, rz, z;
	double *midcolors = (double*)malloc(sizeof(double) * 3);
	double botlshade, topshade, botrshade, lshade, rshade, shade;

	points = OrderPoints(p1, p2, p3, color, shading);

	botl = points[0];
	botlcolor = points[1];
	botlshade = points[6][0];

	top = points[2];
	topcolor = points[3];
	topshade = points[6][1];

	botr = points[4];
	botrcolor = points[5];
	botrshade = points[6][2];

	double rowmin = ceil_441(botl[1]);
	double rowmax = floor_441(top[1]);
	if (rowmax >= 1000)
		rowmax = 999;
	if (rowmin < 0)
		rowmin = 0;

	lslope = slope(botl, top);
	rslope = slope(top, botr);
//printf("ROWMIN: %f, ROWMAX: %f\n", rowmin, rowmax);

	for (double j = rowmin; j <= rowmax; j++)
	{
		if (lslope == 0.0)
			lend = botl[0];
		else
		{
			lb = (top[1] - (lslope*top[0]));
			lend = (j - lb) / lslope;
		}
		if (lend <= 0.0)
			lend = 0.0;
		if (lend > 1000)
			lend = 1000;
		if (rslope == 0.0)
			rend = botr[0];
		else
		{
			rb = (top[1] - (rslope*top[0]));
			rend = (j - rb) / rslope;
		}
		if (rend > 1000)
			rend = 1000;
		if (rend < 0)
			rend = 0;

		lcolors = interpolateColors(botl, botlcolor, top, topcolor, j);
		rcolors = interpolateColors(top, topcolor, botr, botrcolor, j);

		lt = (j - botl[1]) / (top[1] - botl[1]);

		lz = botl[2] + (lt*(top[2] - botl[2]));
		lshade = botlshade + (lt * (topshade - botlshade));

		rt = (j - top[1]) / (botr[1] - top[1]);

		rz = top[2] + (rt*(botr[2] - top[2]));
		rshade = topshade + (rt * (botrshade - topshade));

		for (double k = ceil_441(lend); k <= floor_441(rend); k++)
		{
			t = (k - lend) / (rend - lend);
			midcolors[0] = lcolors[0] + (t*(rcolors[0] - lcolors[0]));
			midcolors[1] = lcolors[1] + (t*(rcolors[1] - lcolors[1]));
			midcolors[2] = lcolors[2] + (t*(rcolors[2] - lcolors[2]));

			z = lz + (t * (rz - lz));
			shade = lshade + (t *(rshade - lshade));

			colorPixel(buffer, midcolors, width, k, j, z, zbuffer, shade);
		}
		free(lcolors);
		free(rcolors);
	}
	free(midcolors);
	free(points);
}

//Find the point that intersects with the top and bottom
void RasterizeArbitrary(double *p1, double *p2, double *p3, unsigned char *buffer, double *zbuffer, double width, double color[3][3], double shading[3])
{
	double y1, y2, y3, m, x, b, t, z;
	int bt, tp, s;
	double *sp, *top, *bottom;
	double newpt[3];
	double newcolors[3][3];
	double newshading[3];

	y1 = p1[1];
	y2 = p2[1];
	y3 = p3[1];
	//Comparisons of y values to see which is the split point, keep track of order because of color change
	if (y1 <= y2 && y2 < y3)
	{
		bottom = p1;
		sp = p2;
		top = p3;
		bt = 0, s = 1, tp = 2;
	}
	else if (y1 <= y3 && y3 < y2)
	{
		bottom = p1;
		sp = p3;
		top = p2;
		bt = 0, s = 2, tp = 1;
	}
	else if (y2 <  y1 && y1 < y3)
	{
		bottom = p2;
		sp = p1;
		top = p3;
		bt = 1, s = 0, tp = 2;
	}
	else if (y2 <= y3 && y3 < y1)
	{
		bottom = p2;
		sp = p3;
		top = p1;
		bt = 1, s = 2, tp = 0;

	}
	else if (y3 < y1 && y1 < y2)
	{
		bottom = p3;
		sp = p1;
		top = p2;
		bt = 2, s = 0, tp = 1;
	}
	else
	{
		bottom = p3;
		sp = p2;
		top = p1;
		bt = 2, s = 1, tp = 0;
	}

	if (bottom[0] > top[0])
		m = slope(top, bottom);
	else
		m = slope(bottom, top);
	if (bottom[0] == top[0])
		x = bottom[0];
	else
	{
		b = (top[1] - (m*top[0]));
		x = (sp[1] - b) / m;
	}
	newpt[0] = x;
	newpt[1] = sp[1];
	t = (sp[1] - bottom[1]) / (top[1] - bottom[1]);		//Interpolate point along known y axis
	z = bottom[2] + (t*(top[2] - bottom[2]));

	newpt[2] = z;

	//fill out correct split point colors
	newcolors[0][0] = color[s][0];
	newcolors[0][1] = color[s][1];
	newcolors[0][2] = color[s][2];
	//calculate correct new point colors
	newcolors[1][0] = color[bt][0] + (t * (color[tp][0] - color[bt][0]));
	newcolors[1][1] = color[bt][1] + (t * (color[tp][1] - color[bt][1]));
	newcolors[1][2] = color[bt][2] + (t * (color[tp][2] - color[bt][2]));
	//For going down, set colors as bottom colors
	newcolors[2][0] = color[bt][0];
	newcolors[2][1] = color[bt][1];
	newcolors[2][2] = color[bt][2];

	newshading[0] = shading[s];
	newshading[1] = shading[bt] + (t*(shading[tp] - shading[bt]));
	newshading[2] = shading[bt];
	
	RasterizeGoingDown(sp, newpt, bottom, buffer, zbuffer, width, newcolors, newshading);
	//For going up, set colors as top colors
	newcolors[2][0] = color[tp][0];
	newcolors[2][1] = color[tp][1];
	newcolors[2][2] = color[tp][2];

	newshading[2] = shading[tp];

	RasterizeGoingUp(sp, newpt, top, buffer, zbuffer, width, newcolors, newshading);
}

double *transformPts(double *pts, Matrix m1, Matrix m2, Matrix m3)
{
	Matrix intermediate;
	Matrix fin;
	double *newpts = (double *)malloc(sizeof(double) * 4);

	//Calculate the final matrix transformation from the 3 matrices
	intermediate = m1.ComposeMatrices(m1, m2);
	fin = m1.ComposeMatrices(intermediate, m3);

	//Transform the old point
	fin.TransformPoint(pts, newpts);

	//Check the q to normalize the x,y,z values
	if (newpts[3] != 0)
	{
		newpts[0] = newpts[0] / newpts[3];
		newpts[1] = newpts[1] / newpts[3];
		newpts[2] = newpts[2] / newpts[3];
		newpts[3] = 0;
	}
	return newpts;
}

double CalculatePhongShading(LightingParameters &, double *viewDirection, double *normal)
{
	double R[3];
	double *L = normalize(lp.lightDir);
	R[0] = 2 * (dotProduct(normal, L)*normal[0]) - L[0];
	R[1] = 2 * (dotProduct(normal, L)*normal[1]) - L[1];
	R[2] = 2 * (dotProduct(normal, L)*normal[2]) - L[2];
	//printf("R = (%f,%f,%f)\n", R[0], R[1], R[2]);

	double *normv = normalize(viewDirection);
	double *normr = normalize(R);
	double btwn = dotProduct(normv, normr);
	//printf("cosine between R and ViewDir is %f\n", btwn);

	double diffuse = (lp.Kd*abs(dotProduct(normal, L)));
	//printf("diffuse = %f\n", diffuse);

	double specular = std::fmax(0, lp.Ks*(pow(btwn, lp.alpha)));
	//printf("specular = %f\n\n\n", specular);


	free(viewDirection);
	free(L);
	free(normv);
	free(normr);
	return (lp.Ka + diffuse + specular);
}

int main()
{
   vtkImageData *image = NewImage(1000, 1000);
   unsigned char *buffer = 
     (unsigned char *) image->GetScalarPointer(0,0,0);
   int npixels = 1000*1000;
   double *zbuffer = (double *)malloc(sizeof(double) * npixels);;
   double p1[4], p2[4], p3[4];
   double shading[3];
   double *pt1, *pt2, *pt3;
   char *typ;
   char filename[32];
   
   for (int k = 0; k < 1; k++)
   {
	   for (int i = 0; i < npixels * 3; i++)
		   buffer[i] = 0;

	   for (int j = 0; j < npixels; ++j)
		   zbuffer[j] = -1;

	   std::vector<Triangle> triangles = GetTriangles();

	   Screen screen;
	   screen.buffer = buffer;
	   screen.width = 1000;
	   screen.height = 1000;

	   Camera c = GetCamera(k, 1000);
	   Matrix m1 = c.CameraTransform();
	   Matrix m2 = c.ViewTransform();
	   Matrix m3 = c.DeviceTransform(1000, 1000);

	   for (int i = 0; i < triangles.size(); ++i)
	   {
		   //Assign triangle points
		   Triangle t = triangles[i];
		   p1[0] = t.X[0];
		   p2[0] = t.X[1];
		   p3[0] = t.X[2];

		   p1[1] = t.Y[0];
		   p2[1] = t.Y[1];
		   p3[1] = t.Y[2];

		   p1[2] = t.Z[0];
		   p2[2] = t.Z[1];
		   p3[2] = t.Z[2];

		   //Add w axis points because we are doing 4D calculations
		   p1[3] = 1;
		   p2[3] = 1;
		   p3[3] = 1;

		   //printf("workind on vertex: (%f, %f, %f)\n", p1[0], p1[1], p1[2]);
		   //printf("normal is: (%f, %f, %f)\n", t.normals[0][0], t.normals[0][1], t.normals[0][2]);
		   shading[0] = CalculatePhongShading(lp, subtract(c.position, p1), t.normals[0]);

		   //printf("workind on vertex: (%f, %f, %f)\n", p2[0], p2[1], p2[2]);
		   //printf("normal is: (%f, %f, %f)\n", t.normals[1][0], t.normals[1][1], t.normals[1][2]);
		   shading[1] = CalculatePhongShading(lp, subtract(c.position, p2), t.normals[1]);

		   //printf("workind on vertex: (%f, %f, %f)\n", p3[0], p3[1], p3[2]);
		   //printf("normal is: (%f, %f, %f)\n", t.normals[2][0], t.normals[2][1], t.normals[2][2]);
		   shading[2] = CalculatePhongShading(lp, subtract(c.position, p3), t.normals[2]);


//printf("Before Transform p1: (%f,%f,%f), p2: (%f,%f,%f), p3: (%f,%f,%f)\n", p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], p3[0], p3[1], p3[2]);
		   //Transform the points according to our new matrix
		   pt1 = transformPts(p1, m1, m2, m3);
		   pt2 = transformPts(p2, m1, m2, m3);
		   pt3 = transformPts(p3, m1, m2, m3);
//printf("After Transform p1: (%f,%f,%f), p2: (%f,%f,%f), p3: (%f,%f,%f)\n", pt1[0], pt1[1], pt1[2], pt2[0], pt2[1], pt2[2], pt3[0], pt3[1], pt3[2]);
		   
		   typ = TriType(pt1, pt2, pt3);

//printf("Triangle: type: %s, p1 = (%f,%f), p2 = (%f,%f), p3 = (%f,%f)\n", typ, p1[0], p1[1], p2[0], p2[1], p3[0], p3[1]);
//printf("Triangle ID: %d\n", i);

		   if (strcmp(typ, "goingdown") == 0)
			   RasterizeGoingDown(pt1, pt2, pt3, buffer, zbuffer, 1000, t.colors, shading);
		   else if (strcmp(typ, "goingup") == 0)
			   RasterizeGoingUp(pt1, pt2, pt3, buffer, zbuffer, 1000, t.colors, shading);
		   else
			   RasterizeArbitrary(pt1, pt2, pt3, buffer, zbuffer, 1000, t.colors, shading);

		   free(pt1);
		   free(pt2);
		   free(pt3);
	   }

	   //check which camera you are dealing with and name accordingly
	   sprintf(filename, "frame%03d", k);
	   WriteImage(image, filename);
	   printf("Wrote image %s\n", filename);
   }
   free(buffer);
   free(zbuffer);
}