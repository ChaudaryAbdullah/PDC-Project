#ifndef OPENCL_UTILS_H
#define OPENCL_UTILS_H

#include <CL/cl.h>
#include <vector>
#include <string>

struct OpenCLContext
{
    cl_platform_id platform;
    std::vector<cl_device_id> devices;
    cl_context context;
    cl_program program;
    cl_command_queue queue;
};

bool setupOpenCL(OpenCLContext &ctx, const std::string &kernel_file);
void cleanupOpenCL(OpenCLContext &ctx);
bool runRelaxationKernel(OpenCLContext &ctx,
                         std::vector<float> &dist,
                         std::vector<int> &parent,
                         const std::vector<std::pair<int, int>> &edges,
                         const std::vector<float> &weights);

#endif // OPENCL_UTILS_H