import subprocess
import time
import matplotlib.pyplot as plt
import os
import re
import glob
import math

def extract_number(filename):
    # Extract the number from the filename (e.g., "update1.txt" -> 1)
    match = re.search(r'update(\d+)\.txt', filename)
    return int(match.group(1)) if match else float('inf')

def run_serial_program(input_file, update_file, source_vertex, output_file):
    if not os.path.exists(input_file) or not os.path.exists(update_file):
        print(f"Error: Input file {input_file} or update file {update_file} does not exist.")
        return None, 0

    print(f"Running serial SSSP with update file {update_file}...")

    start_time = time.time()
    
    command = ['./serial_sssp', input_file, update_file, str(source_vertex), output_file]
    print(f"Executing serial command: {' '.join(command)}")
    
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()
    
    if process.returncode != 0:
        print(f"Error running serial SSSP with {update_file}: {error.decode()}")
        print(f"Serial standard output: {output.decode()}")
        return None, 0
    
    end_time = time.time()
    execution_time = (end_time - start_time) * 1000  # Convert to milliseconds
    
    output_str = output.decode()
    match = re.search(r"SSSP update completed in (\d+) seconds", output_str)
    if match:
        execution_time = int(match.group(1)) * 1000  # Convert to milliseconds
        print(f"Serial SSSP reported execution time: {execution_time} ms")
    else:
        print(f"Serial SSSP execution time (measured): {execution_time:.2f} ms")
    
    # Estimate number of updates by counting lines (approximate)
    update_count = sum(1 for line in open(update_file) if line.strip() and not line.startswith('#'))
    return execution_time, update_count

def run_parallel_program(input_file, update_file, source_vertex, num_processes, output_file, use_openmp=True, use_opencl=True):
    if not os.path.exists(input_file) or not os.path.exists(update_file):
        print(f"Error: Input file {input_file} or update file {update_file} does not exist.")
        return None, 0

    print(f"Running parallel SSSP with {num_processes} processes and update file {update_file}...")

    start_time = time.time()
    
    command = ['mpirun', '--use-hwthread-cpus', '--bind-to', 'core:overload-allowed', '-np', str(num_processes), './sssp', input_file, update_file, str(source_vertex), output_file]
    if use_openmp:
        command.append('--openmp')
    if use_opencl:
        command.append('--opencl')
    
    print(f"Executing parallel command: {' '.join(command)}")
    
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()
    
    if process.returncode != 0:
        print(f"Error running parallel SSSP with {num_processes} processes and {update_file}: {error.decode()}")
        print(f"Parallel standard output: {output.decode()}")
        return None, 0
    
    end_time = time.time()
    execution_time = (end_time - start_time) * 1000  # Convert to milliseconds
    print(f"Completed parallel run with {num_processes} processes in {execution_time:.2f} ms.")
    
    # Estimate number of updates by counting lines (approximate)
    update_count = sum(1 for line in open(update_file) if line.strip() and not line.startswith('#'))
    return execution_time, update_count

def compare_serial_parallel(input_file, update_file, source_vertex, num_processes_list=[1, 2, 4, 6, 8], use_openmp=True, use_opencl=True):
    serial_output_file = f"output_serial_{os.path.basename(update_file).replace('.txt', '')}.txt"
    serial_time, serial_updates = run_serial_program(input_file, update_file, source_vertex, serial_output_file)
    
    parallel_results = []
    for num_processes in num_processes_list:
        parallel_output_file = f"output_parallel_{os.path.basename(update_file).replace('.txt', '')}_{num_processes}.txt"
        execution_time, _ = run_parallel_program(input_file, update_file, source_vertex, num_processes, parallel_output_file, use_openmp, use_opencl)
        if execution_time is not None:
            parallel_results.append((num_processes, execution_time))
        else:
            print(f"Skipping parallel run with {num_processes} processes due to execution failure.")
    
    return serial_time, serial_updates, parallel_results

def plot_comparison(results_dict):
    if not results_dict:
        print("No results to plot due to execution failures.")
        return

    # Prepare data for plotting
    update_counts = []
    serial_times = []
    parallel_times = []
    for update_file, (serial_time, update_count, parallel_results) in results_dict.items():
        if serial_time is not None and parallel_results:
            update_counts.append(update_count)
            serial_times.append(serial_time)
            # Take the first parallel result (for num_processes=4)
            parallel_times.append(parallel_results[0][1])  # [1] is execution_time

    if not update_counts or not serial_times or not parallel_times:
        print("Insufficient data for plotting.")
        return

    # Sort data by update_counts to ensure monotonic x-axis
    sorted_data = sorted(zip(update_counts, serial_times, parallel_times), key=lambda x: x[0])
    update_counts, serial_times, parallel_times = zip(*sorted_data)

    # Debug: Print data to verify
    print("Plotting data:")
    for uc, st, pt in zip(update_counts, serial_times, parallel_times):
        print(f"Updates: {uc}, Serial Time: {st} ms, Parallel Time: {pt} ms")

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), gridspec_kw={'height_ratios': [1, 1]})

    # Bar Chart (Top Plot)
    x = range(len(update_counts))
    width = 0.35
    ax1.bar([i - width/2 for i in x], serial_times, width, label='Sequential (ms)', color='blue')
    ax1.bar([i + width/2 for i in x], parallel_times, width, label='Hybrid (ms)', color='orange')
    ax1.set_ylabel('Time (ms)')
    ax1.set_title('Comparison of Sequential vs Hybrid Processing Time')
    ax1.set_xticks(x)
    ax1.set_xticklabels(update_counts)
    ax1.legend()
    ax1.grid(True)

    # Line Chart (Bottom Plot)
    ax2.plot(update_counts, serial_times, marker='o', label='Sequential (ms)', color='blue')
    ax2.plot(update_counts, parallel_times, marker='o', label='Hybrid (ms)', color='orange')
    ax2.set_xscale('linear')  # Use linear scale for linear appearance
    ax2.set_xlabel('Number of Updates')
    ax2.set_ylabel('Time (ms)')
    ax2.set_title('Performance Scaling with Increasing Updates')
    ax2.legend()
    ax2.grid(True)

    plt.tight_layout()

    try:
        plt.savefig('sssp_comparison.png')
        print("Comparison plot saved as sssp_comparison.png")
    except Exception as e:
        print(f"Failed to save plot: {e}")
        plt.show()

    plt.close()

def main():
    input_file = "sample_graph.txt"
    source_vertex = 10000
    
    update_files = glob.glob("sample_update*.txt")
    if not update_files:
        print("No update files found matching 'sample_update*.txt'.")
        return
    
    # Sort update files by numerical suffix
    update_files.sort(key=extract_number)
    print(f"Found {len(update_files)} update files: {update_files}")
    
    results_dict = {}
    for update_file in update_files:
        print(f"\nProcessing update file: {update_file}")
        serial_time, update_count, parallel_results = compare_serial_parallel(
            input_file, update_file, source_vertex, 
            num_processes_list=[4], 
            use_openmp=True, use_opencl=True
        )
        results_dict[update_file] = (serial_time, update_count, parallel_results)
    
    plot_comparison(results_dict)

if __name__ == "__main__":
    main()