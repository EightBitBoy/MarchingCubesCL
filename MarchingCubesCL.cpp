#include <stdexcept>
#include <utility>
#include <fstream>
#include <fantom/algorithm.hpp>
#include <fantom/graphics.hpp>
#include <fantom/fields.hpp>
#include <fantom/gui.hpp>
#include <CL/cl.hpp>
#include "MarchingCubesKernel.hpp"
#include "MarchingCubesTables.hpp"

#define __CL_ENABLE_EXCEPTIONS

using namespace std;
using namespace fantom;

namespace MC
{
	class MarchingCubes: public Algorithm
	{
	public:
		unique_ptr<Graphics> polygonGroup;

		struct Options: public Algorithm::Options
		{
			Options()
			{
				add<const TensorField<3, Scalar>*>("field", "", 0);
				add<const Grid<3>*>("grid", "", 0);
				add<float>("isoValue", "the iso value used by the Marching Cubes algorithm", 0.00080);
				add<int>("numberOfUsedCells", "", 100);
				add<Color>("color", "", Color(0.75, 0.0, 0.0));
			}
		};

		struct MyWindow
		{

		  DockWindow mWindow;
		  BoxLayout mLayout;
		  LineEdit mText;
		  PushButton mButton;
		  MarchingCubes& mAlgo;

		  MyWindow(MainWindow& mainWindow, MarchingCubes& algo)
			: mWindow(mainWindow, DockWindow::FFREE, "algorithm window"),
			  mLayout(mWindow.getWidgetHolder(), false),
			  mText(mLayout.addWidgetHolder()),
			  mButton(mLayout.addWidgetHolder(), "Send", bind(&MyWindow::send, this)),
			  mAlgo(algo)
		  {
		  }

		  void send()
		  {
			//mAlgo.scheduleJob(bind(&MarchingCubes::printError, &mAlgo, mText.get()));
		  }

		};

		void printError(cl_int error, string s)
		{
			if(error != CL_SUCCESS)
			{
				debugLog() << "ERROR: " << error << " " << s << endl;
			}
		}

		Point3 interpolate(Point3 pointA, Point3 pointB, float valueA, float valueB, float isoValue)
		{
			Point3 result;

			float factor;
			factor = (isoValue - valueA) / (valueB - valueA);
			result[0] = pointA[0] + factor * (pointB[0] - pointA[0]);
			result[1] = pointA[1] + factor * (pointB[1] - pointA[1]);
			result[2] = pointA[2] + factor * (pointB[2] - pointA[2]);
			return result;
		}

		MarchingCubes(const Parameters& parameters): Algorithm(parameters)
		{
			// foo
			// get the data
			const TensorField<3, Scalar>* field = parameters.get<const TensorField<3, Scalar>*>("field");
			const Grid<3>* grid = parameters.get<const Grid<3>*>("grid");
			const float isoValue = parameters.get<float>("isoValue");
			const int numberOfUsedCells = parameters.get<int>("numberOfUsedCells");
			const Color color = parameters.get<Color>("color");

			if(field == false)
			{
				throw runtime_error("field not set!");
			}
			if(grid == false)
			{
				throw runtime_error("grid not set!");
			}

			// add some test objects
			polygonGroup = makeGraphics("polygonGroup");
			//polygonGroup->primitive().addSphere(Point3(0, 0, 0), 0.5, color);

			// some useful variables
			cl_int error = CL_SUCCESS;
			cl_ulong infoNumber = 0;
			string infoString = "";

			// Create the two input vectors
			const int LIST_SIZE = 10;
			int *A = new int[LIST_SIZE]; 
			int *B = new int[LIST_SIZE];
			for(int i = 0; i < LIST_SIZE; i++) {
				A[i] = i;
				B[i] = LIST_SIZE - i;
			}
 
			try { 
				// get the platforms
				vector<cl::Platform> platforms;
				cl::Platform::get(&platforms);
				for(int i = 0; i < platforms.size(); i++)
				{
					platforms[i].getInfo(CL_PLATFORM_NAME, &infoString);
					debugLog() << "platform #" << i << " name: " << infoString << endl;
					platforms[i].getInfo(CL_PLATFORM_VERSION, &infoString);
					debugLog() << "platform #" << i << " version: " << infoString << endl;
				}
				debugLog() << "" << endl;

				// select the default platform 0, create a context, use the GPU
				cl_context_properties properties[3] = { 
					CL_CONTEXT_PLATFORM, 
					(cl_context_properties)(platforms[0])(), 
					0 
				};
				cl::Context context(CL_DEVICE_TYPE_GPU, properties, NULL, NULL, &error);
				printError(error, "context");
 
				// get the devices
				vector<cl::Device> devices = context.getInfo<CL_CONTEXT_DEVICES>();
				if(devices.size() == 0)
				{
					throw runtime_error("No devices found!");
				}
				for(int i = 0; i < devices.size(); i++)
				{
					devices[i].getInfo(CL_DEVICE_NAME, &infoString);
					debugLog() << "device #"<< i << " name: " << infoString << endl;
					devices[i].getInfo(CL_DEVICE_GLOBAL_MEM_SIZE, &infoNumber);
					debugLog() << "device #"<< i << " global memory size (MB): " << (infoNumber/1024/1024) << endl;
					devices[i].getInfo(CL_DEVICE_GLOBAL_MEM_CACHE_SIZE, &infoNumber);
					debugLog() << "device #"<< i << " global memory cache size (MB): " << (infoNumber/1024/1024) << endl;
					devices[i].getInfo(CL_DEVICE_MAX_CLOCK_FREQUENCY, &infoNumber);
					debugLog() << "device #"<< i << " max clock frequency: " << infoNumber << endl;
					devices[i].getInfo(CL_DEVICE_MAX_COMPUTE_UNITS, &infoNumber);
					debugLog() << "device #"<< i << " max compute units: " << infoNumber << endl;
					devices[i].getInfo(CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, &infoNumber);
					debugLog() << "device #"<< i << " max constant buffer size (KB): " << (infoNumber/1024) << endl;
				}
				debugLog() << "" << endl;
				// create the command queue, use the first device
				cl::CommandQueue queue = cl::CommandQueue(context, devices[0]);
 
				ifstream file("code.cl");
				string prog(istreambuf_iterator<char>(file), (istreambuf_iterator<char>()));





				// get the kernel source code
				string s = getKernelSource();
				cl::Program::Sources source(1, std::make_pair(prog.c_str(), prog.length()+1));

				// Make program of the source code in the context
				cl::Program program = cl::Program(context, source);
 
				// Build program for these specific devices
				error = program.build(devices);
				string f;
				program.getBuildInfo(devices[0], CL_PROGRAM_BUILD_LOG, &f);
				debugLog() << f << endl;
				printError(error, "kernel");
				// Make kernel
				cl::Kernel kernel(program, "vector_add");


				////////////////////////////////////////
				////////////////////////////////////////

				// prepare the data
				auto evaluator = field->makeEvaluator();
				auto& points = grid->parent().points();

				debugLog() << "cells: " << grid->numCells() << endl;
				//const size_t numberOfCells = grid->numCells(); // TODO use the correct number
				//const size_t numberOfCells = 100;
				const size_t numberOfCells = numberOfUsedCells;
				const int numberOfCellPoints = 8;
				const size_t numberOfValues = numberOfCells * numberOfCellPoints;
				float *values = new float[numberOfCells];

				//float values[numberOfValues];
				double time = 0;

				for(Progress i(*this, "working, size", numberOfCells); i < numberOfCells; ++i)
				{
					Cell cell = grid->cell(i);
					for(size_t j = 0; j < numberOfCellPoints; ++j)
					{
						Point3 p = points[cell.index(j)];
						if(evaluator->reset(points[cell.index(j)], time))
						{
							values[(i * numberOfCellPoints) + j] = evaluator->value()();
						}
						else
						{
							values[(i * numberOfCellPoints) + j] = 0;
						}
					}
				}
				
				////////////////////////////////////////
				////////////////////////////////////////


				float valuesTest[numberOfValues];
				for(int i = 0; i < numberOfValues; i++)
				{
					valuesTest[i] = 99.8877;
				}

				// Create memory buffers
				cl::Buffer bufferA = cl::Buffer(context, CL_MEM_READ_ONLY, LIST_SIZE * sizeof(int));
				cl::Buffer bufferB = cl::Buffer(context, CL_MEM_READ_ONLY, LIST_SIZE * sizeof(int));
				cl::Buffer bufferC = cl::Buffer(context, CL_MEM_WRITE_ONLY, LIST_SIZE * sizeof(int));

				cl::Buffer bufferIsoValue = cl::Buffer(context, CL_MEM_READ_ONLY, 1 * sizeof(float));
				cl::Buffer bufferEdgeList = cl::Buffer(context, CL_MEM_READ_ONLY, 256 * sizeof(int));
				cl::Buffer bufferValues = cl::Buffer(context, CL_MEM_READ_ONLY, numberOfValues * sizeof(float));
				cl::Buffer bufferOutputInterpolatedValues = cl::Buffer(context, CL_MEM_WRITE_ONLY, numberOfCells * 12 * sizeof(float));
				cl::Buffer bufferOutputCubeIndices = cl::Buffer(context, CL_MEM_WRITE_ONLY, numberOfCells * sizeof(int));
 
				// Copy lists A and B to the memory buffers
				queue.enqueueWriteBuffer(bufferA, CL_TRUE, 0, LIST_SIZE * sizeof(int), A);
				queue.enqueueWriteBuffer(bufferB, CL_TRUE, 0, LIST_SIZE * sizeof(int), B);

				queue.enqueueWriteBuffer(bufferIsoValue, CL_TRUE, 0, 1 * sizeof(float), &isoValue);
				queue.enqueueWriteBuffer(bufferValues, CL_TRUE, 0, numberOfValues * sizeof(float), valuesTest);
				queue.enqueueWriteBuffer(bufferEdgeList, CL_TRUE, 0, 256 * sizeof(int), edgeTable);
 

				// Set arguments to kernel
				kernel.setArg(0, bufferA);
				kernel.setArg(1, bufferB);
				kernel.setArg(2, bufferC);

				error = kernel.setArg(3, isoValue);
				printError(error, "isoValue");
				kernel.setArg(4, bufferEdgeList);
				kernel.setArg(5, bufferValues);
				kernel.setArg(6, bufferOutputInterpolatedValues);
				kernel.setArg(7, bufferOutputCubeIndices);
 
				// Run the kernel on specific ND range
				cl::NDRange global(LIST_SIZE); // TODO use the number of cells
				cl::NDRange local(1);
				queue.enqueueNDRangeKernel(kernel, cl::NullRange, global, local);
 
				// Read buffer C into a local list

				int C[LIST_SIZE];
				queue.enqueueReadBuffer(bufferC, CL_TRUE, 0, LIST_SIZE * sizeof(int), C);
				for(int i = 0; i < LIST_SIZE; i++)
				{
					//debugLog() << A[i] << " + " << B[i] << " = " << C[i] << endl;
				}

				float interpolatedValues[numberOfCells * 12];
				queue.enqueueReadBuffer(bufferOutputInterpolatedValues, CL_TRUE, 0, numberOfCells * 12 * sizeof(float), interpolatedValues);
				for(int i = 0; i < numberOfValues; i++)
				{
					//debugLog() << "out: " << interpolatedValues[i] << endl;
				}

				int cubeIndices[numberOfCells];
				queue.enqueueReadBuffer(bufferOutputCubeIndices, CL_TRUE, 0, numberOfCells * sizeof(int), cubeIndices);
				for(int i = 0; i < numberOfCells; i++)
				{
					//debugLog() << cubeIndices[i] << endl;
				}

				
				// render the triangles
				for(Progress i(*this, "working, size", numberOfCells); i < numberOfCells; ++i)
				{
					debugLog() << "working" << endl;
					Cell cell = grid->cell(i);

					Point3 cubePoints[8];
					for(size_t j = 0; j < numberOfCellPoints; ++j)
					{
						cubePoints[j] = points[cell.index(j)];
					}

					int index = cubeIndices[i];
					Point3 vertices[12];
					if(edgeTable[index] == 0)
						continue;
					

					if(edgeTable[index] & 1)
						vertices[0] = interpolate(cubePoints[0], cubePoints[1], values[(8 * i) + 0], values[(8 * i) + 1], isoValue);
					if(edgeTable[index] & 2)
						vertices[1] = interpolate(cubePoints[1], cubePoints[2], values[(8 * i) + 1], values[(8 * i) + 2], isoValue);
					if(edgeTable[index] & 4)
						vertices[2] = interpolate(cubePoints[2], cubePoints[3], values[(8 * i) + 2], values[(8 * i) + 3], isoValue);
					if(edgeTable[index] & 8)
						vertices[3] = interpolate(cubePoints[3], cubePoints[0], values[(8 * i) + 3], values[(8 * i) + 0], isoValue);
					if(edgeTable[index] & 16)
						vertices[4] = interpolate(cubePoints[4], cubePoints[5], values[(8 * i) + 4], values[(8 * i) + 5], isoValue);
					if(edgeTable[index] & 32)
						vertices[5] = interpolate(cubePoints[5], cubePoints[6], values[(8 * i) + 5], values[(8 * i) + 6], isoValue);
					if(edgeTable[index] & 64)
						vertices[6] = interpolate(cubePoints[6], cubePoints[7], values[(8 * i) + 6], values[(8 * i) + 7], isoValue);
					if(edgeTable[index] & 128)
						vertices[7] = interpolate(cubePoints[7], cubePoints[4], values[(8 * i) + 7], values[(8 * i) + 4], isoValue);
					if(edgeTable[index] & 256)
						vertices[8] = interpolate(cubePoints[0], cubePoints[4], values[(8 * i) + 0], values[(8 * i) + 4], isoValue);
					if(edgeTable[index] & 512)
						vertices[9] = interpolate(cubePoints[1], cubePoints[5], values[(8 * i) + 1], values[(8 * i) + 5], isoValue);
					if(edgeTable[index] & 1024)
						vertices[10] = interpolate(cubePoints[2], cubePoints[6], values[(8 * i) + 2], values[(8 * i) + 6], isoValue);
					if(edgeTable[index] & 2048)
						vertices[11] = interpolate(cubePoints[3], cubePoints[7], values[(8 * i) + 3], values[(8 * i) + 7], isoValue);

					for(int j = 0; triTable[index][j] != -1; j += 3)
					{
						debugLog() << "YEAH" << endl;
						vector<Point3> polygonPoints;
						polygonPoints.push_back(vertices[j + 0]);
						polygonPoints.push_back(vertices[j + 1]);
						polygonPoints.push_back(vertices[j + 2]);

						polygonGroup->primitive().add(Primitive::POLYGON).setColor(color).setVertices(polygonPoints);
					}
				}
			}
			catch(...)
			{
				throw runtime_error("SOMETHING WENT WRONG!!!");
			}
		}
	};

	AlgorithmRegister<MarchingCubes> dummy("MarchingCubesCL", "executes MarchingCubesCL");
}
