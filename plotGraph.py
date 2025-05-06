import subprocess
import time
import matplotlib.pyplot as plt
import os

def run_cpp_program(input_file, update_file, num_processes=1, use_openmp=False, use_opencl=False):
    # Check if input files exist
    if not os.path.exists(input_file):
        print(f"Error: Input file {input_file} does not exist.")
        return None
    if not os.path.exists(update_file):
        print(f"Error: Update file {update_file} does not exist.")
        return None

    print(f"Running sssp with {num_processes} processes...")

    start_time = time.time()
    
    # Command to run the MPI program with --use-hwthread-cpus and --bind-to core:overload-allowed
    command = ['mpirun', '--use-hwthread-cpus', '--bind-to', 'core:overload-allowed', '-np', str(num_processes), './sssp', input_file, update_file, '10000', 'output.txt']
    if use_openmp:
        command.append('--openmp')
    if use_opencl:
        command.append('--opencl')
    
    print(f"Executing command: {' '.join(command)}")
    
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()
    
    if process.returncode != 0:
        print(f"Error running {input_file} with {num_processes} processes: {error.decode()}")
        print(f"Standard output: {output.decode()}")
        return None
    
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Completed run with {num_processes} processes in {execution_time:.2f} seconds.")
    
    return execution_time

def test_multiple_processes(input_file, update_file, num_processes_list=[1, 2, 4, 6, 8], use_openmp=False, use_opencl=False):
    results = []
    
    for num_processes in num_processes_list:
        execution_time = run_cpp_program(input_file, update_file, num_processes, use_openmp, use_opencl)
        if execution_time is not None:
            results.append((num_processes, execution_time))
        else:
            print(f"Skipping {num_processes} processes due to execution failure.")
    
    return results

def plot_results(results, type='bar'):
    if not results:
        print("No results to plot due to execution failures.")
        return
    
    num_processes_list = [r[0] for r in results]
    execution_times = [r[1] for r in results]
    
    plt.figure(figsize=(8, 6))
    if type == 'bar':
        plt.bar([str(np) for np in num_processes_list], execution_times, color='skyblue', alpha=0.8)
        plt.xlabel('Number of MPI Processes')
        plt.ylabel('Execution Time (seconds)')
        plt.title('Execution Time of SSSP on Sample Graph with Different MPI Processes')
    elif type == 'line':
        plt.plot(num_processes_list, execution_times, marker='o', color='skyblue')
        plt.xlabel('Number of MPI Processes')
        plt.ylabel('Execution Time (seconds)')
        plt.title('Execution Time of SSSP on Sample Graph with Different MPI Processes')
    elif type == 'scatter':
        plt.scatter(num_processes_list, execution_times, color='skyblue')
        plt.xlabel('Number of MPI Processes')
        plt.ylabel('Execution Time (seconds)')
        plt.title('Execution Time of SSSP on Sample Graph with Different MPI Processes')
    else:
        raise ValueError("Unsupported plot type. Use 'bar', 'line', or 'scatter'.")
    
    plt.tight_layout()
    
    try:
        plt.savefig('sssp_performance.png')
        print("Plot saved as sssp_performance.png")
    except Exception as e:
        print(f"Failed to save plot: {e}")
        plt.show()  # Fallback to display the plot
    
    plt.close()

def main():
    # Fixed input files
    input_file = "sample_graph.txt"
    update_file = "sample_updates.txt"
    
    # Test with different numbers of MPI processes, OpenMP and OpenCL enabled
    results = test_multiple_processes(input_file, update_file, num_processes_list=[2, 3, 4], use_openmp=True, use_opencl=True)
    
    plot_results(results, 'bar')

if __name__ == "__main__":
    main()