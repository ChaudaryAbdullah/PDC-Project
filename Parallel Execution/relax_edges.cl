// Optimized kernel for edge relaxation in SSSP

#define INF INFINITY

__kernel void relax_edges(__global float *dist, __global int *parent,
                          __global int2 *edges, __global float *weights,
                          const int num_edges)
{
    int i = get_global_id(0);
    if (i >= num_edges)
        return;

    // Local variables for better performance
    int u = edges[i].x;
    int v = edges[i].y;
    float w = weights[i];
    float dist_u = dist[u];
    float dist_v = dist[v];

    // Check for infinities
    if (dist_u == INF)
        return;

    // Calculate potential new distance
    float new_dist = dist_u + w;

    // Perform relaxation if better path found
    if (new_dist < dist_v)
    {
        // Use atomic operation for thread safety
        atomic_min_float(&dist[v], new_dist);

        // Check if we actually updated the distance before changing parent
        if (dist[v] == new_dist)
        {
            parent[v] = u;
        }
    }
}

// Custom atomic min operation for floats since OpenCL doesn't provide it natively
inline void atomic_min_float(__global float *addr, float val)
{
    union
    {
        float f;
        unsigned int i;
    } old_val, new_val;

    do
    {
        old_val.f = *addr;
        new_val.f = (val < old_val.f) ? val : old_val.f;
    } while (atomic_cmpxchg((volatile __global unsigned int *)addr,
                            old_val.i, new_val.i) != old_val.i);
}