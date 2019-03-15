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

int counter = 0;

using std::cerr;
using std::endl;


double ceil_1(double f)
{
    return ceil(f-0.00001);
}

double 
_1(double f)
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

class Triangle
{
  public:
      double         X[3];
      double         Y[3];
double		 Z[3];
double		 colors[3][3];

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
	rdr->SetFileName("proj1d_geometry.vtk");
	cerr << "Reading" << endl;
	rdr->Update();
	cerr << "Done reading" << endl;
	if (rdr->GetOutput()->GetNumberOfCells() == 0)
	{
		cerr << "Unable to open file!!" << endl;
		exit(EXIT_FAILURE);
	}
	vtkPolyData *pd = rdr->GetOutput();
	int numTris = pd->GetNumberOfCells();
	vtkPoints *pts = pd->GetPoints();
	vtkCellArray *cells = pd->GetPolys();
	vtkFloatArray *var = (vtkFloatArray *)pd->GetPointData()->GetArray("hardyglobal");
	float *color_ptr = var->GetPointer(0);
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
		tris[idx].X[0] = pts->GetPoint(ptIds[0])[0];
		tris[idx].X[1] = pts->GetPoint(ptIds[1])[0];
		tris[idx].X[2] = pts->GetPoint(ptIds[2])[0];
		tris[idx].Y[0] = pts->GetPoint(ptIds[0])[1];
		tris[idx].Y[1] = pts->GetPoint(ptIds[1])[1];
		tris[idx].Y[2] = pts->GetPoint(ptIds[2])[1];
		tris[idx].Z[0] = pts->GetPoint(ptIds[0])[2];
		tris[idx].Z[1] = pts->GetPoint(ptIds[1])[2];
		tris[idx].Z[2] = pts->GetPoint(ptIds[2])[2];
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
double **OrderPoints(double *p1, double *p2, double *p3, double colors[3][3])
{
	double **points = (double**)malloc(sizeof(double*) * 6);
	for (int i = 0; i < 4; ++i)
	{
		points[i] = (double*)malloc(sizeof(double) * 3);
	}
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
		}
		else
		{
			points[0] = p2;
			points[1] = colors[1];
			points[4] = p1;
			points[5] = colors[0];
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
		}
		else
		{
			points[0] = p3;
			points[1] = colors[2];
			points[4] = p2;
			points[5] = colors[1];
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
	}
	else
	{
		points[0] = p3;
		points[1] = colors[2];
		points[2] = p2;
		points[3] = colors[1];
		points[4] = p1;
		points[5] = colors[0];
	}
	return points;
}

double slope(double *pt1, double *pt2)
{
	if (pt1[0] == pt2[0])
		return 0.0;
	return ((pt2[1] - pt1[1]) / (pt2[0] - pt1[0]));
}


void colorPixel(unsigned char *buffer, double *colors, int width, int x, int y, double z, double *zbuffer)
{
	//check if pixels are trying to overlap onto left side
	if (x == 1000 && y >= 0)
		return;
	//check depth buffer to see if we should display these
	if (z >= zbuffer[(width*y) + x])
	{
		zbuffer[(width*y) + x] = z;
		buffer[(3 * width*y) + 3 * x] = ceil_1(colors[0]*255);
		buffer[(3 * width*y) + 3 * x + 1] = ceil_1(colors[1] * 255);
		buffer[(3 * width*y) + 3 * x + 2] = ceil_1(colors[2] * 255);
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
void RasterizeGoingDown(double *p1, double *p2, double *p3, unsigned char *buffer, double *zbuffer, double width, double color[3][3])
{
//printf("Inside Going Down Rasterize\n");
	double **points;
	double *topl, *topr, *bottom;
	double *toplcolor, *toprcolor, *bottomcolor;
	double lslope, rslope, lend, rend, lb, rb;
	double *lcolors, *rcolors;
	double *midcolors = (double*)malloc(sizeof(double) * 3);
	double lt, rt, t, lz, rz, z;

	//return points in order followed by color
	points = OrderPoints(p1, p2, p3, color);
	topl = points[0];
	toplcolor = points[1];
	bottom = points[2];
	bottomcolor = points[3];
	topr = points[4];
	toprcolor = points[5];

	double rowmin = ceil_1(bottom[1]);
	double rowmax = floor_1(topl[1]);
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

		//Interpolate the z value for the rightend
		rt = (j - bottom[1]) / (topr[1] - bottom[1]);
		rz = bottom[2] + (rt*(topr[2] - bottom[2]));
		
//printf("lz = %f + (%f * (%f - %f))\n", topl[2], lt, bottom[2], topl[2]);
//printf("rz = %f + (%f * (%f - %f))\n", bottom[2], rt, topr[2], bottom[2]);

		//Go through all points along the line
		for (double k = ceil_1(lend); k <= floor_1(rend); k++)
		{
			//calculate colors for the individual pixel
			t = (k - lend) / (rend - lend);
			midcolors[0] = lcolors[0] + (t*(rcolors[0] - lcolors[0]));
			midcolors[1] = lcolors[1] + (t*(rcolors[1] - lcolors[1]));
			midcolors[2] = lcolors[2] + (t*(rcolors[2] - lcolors[2]));

			//calculate depth to pass to colorPixel
			z = lz + (t * (rz - lz));
			//color pixel
			colorPixel(buffer, midcolors, width, k, j, z, zbuffer);
		}
	}
}
//Calculate pixels for any going up triangles
void RasterizeGoingUp(double *p1, double *p2, double *p3, unsigned char *buffer, double *zbuffer, double width, double color[3][3])
{
//printf("Inside Going UP Rasterize\n");
	double **points;
	double *botl, *botr, *top;
	double lslope, rslope, lend, rend, lb, rb, t;
	double *botlcolor, *topcolor, *botrcolor;
	double *lcolors, *rcolors;
	double lt, rt, lz, rz, z;
	double *midcolors = (double*)malloc(sizeof(double) * 3);

	points = OrderPoints(p1, p2, p3, color);

	botl = points[0];
	botlcolor = points[1];
	top = points[2];
	topcolor = points[3];
	botr = points[4];
	botrcolor = points[5];

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

		rt = (j - top[1]) / (botr[1] - top[1]);
		rz = top[2] + (rt*(botr[2] - top[2]));

		for (double k = ceil_441(lend); k <= floor_1(rend); k++)
		{
			t = (k - lend) / (rend - lend);
			midcolors[0] = lcolors[0] + (t*(rcolors[0] - lcolors[0]));
			midcolors[1] = lcolors[1] + (t*(rcolors[1] - lcolors[1]));
			midcolors[2] = lcolors[2] + (t*(rcolors[2] - lcolors[2]));
			z = lz + (t * (rz - lz));

			colorPixel(buffer, midcolors, width, k, j, z, zbuffer);
		}
	}
}

//Find the point that intersects with the top and bottom
void RasterizeArbitrary(double *p1, double *p2, double *p3, unsigned char *buffer, double *zbuffer, double width, double color[3][3])
{
	double y1, y2, y3, m, x, b, t, z;
	int bt, tp, s;
	double *sp, *top, *bottom;
	double newpt[3];
	double newcolors[3][3];

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
	
	RasterizeGoingDown(sp, newpt, bottom, buffer, zbuffer, width, newcolors);
	//For going up, set colors as top colors
	newcolors[2][0] = color[tp][0];
	newcolors[2][1] = color[tp][1];
	newcolors[2][2] = color[tp][2];

	RasterizeGoingUp(sp, newpt, top, buffer, zbuffer, width, newcolors);
}


int main()
{
   vtkImageData *image = NewImage(1000, 1000);
   unsigned char *buffer = 
     (unsigned char *) image->GetScalarPointer(0,0,0);
   int npixels = 1000*1000;
   for (int i = 0 ; i < npixels*3 ; i++)
       buffer[i] = 0;

   double *zbuffer = (double *)malloc(sizeof(double) * npixels);;
   for (int j = 0; j < npixels; ++j)
	   zbuffer[j] = -1;

   std::vector<Triangle> triangles = GetTriangles();

   Screen screen;
   screen.buffer = buffer;
   screen.width = 1000;
   screen.height = 1000;

   double p1[3], p2[3], p3[3];
   char *typ;
   
//printf("NUMTRI = %d\n", triangles.size());
   for (int i = 0; i < triangles.size(); ++i)
   {
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
	   typ = TriType(p1, p2, p3);
	   //printf("Triangle: type: %s, p1 = (%f,%f), p2 = (%f,%f), p3 = (%f,%f)\n", typ, p1[0], p1[1], p2[0], p2[1], p3[0], p3[1]);
	   //printf("Triangle ID: %d\n", i);

	   if (strcmp(typ, "goingdown") == 0)
		   RasterizeGoingDown(p1, p2, p3, buffer, zbuffer, 1000, t.colors);
	   else if (strcmp(typ, "goingup") == 0)
		   RasterizeGoingUp(p1, p2, p3, buffer, zbuffer, 1000, t.colors);
	   else
		   RasterizeArbitrary(p1, p2, p3, buffer, zbuffer, 1000, t.colors);
   }
   WriteImage(image, "allTriangles");
}
