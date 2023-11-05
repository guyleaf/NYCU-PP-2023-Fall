# Programming Assignment II: Multi-thread Programming
**Parallel Programming by Prof. Yi-Ping You**

Assignment description: https://pp-f23.github.io/assignments/HW2

## Introduction
The purpose of this assignment is to familiarize yourself with Pthread and std::thread programming in C and C++, respectively. You will also gain experience measuring and reasoning about the performance of parallel programs (a challenging, but important, skill you will use throughout this class). This assignment involves only a small amount of programming, but a lot of analysis!

## Q1
### Question
In your write-up, produce a graph of **speedup compared to the reference sequential implementation** as a function of the number of threads used **FOR VIEW 1**. 

Is speedup linear in the number of threads used?  
In your writeup hypothesize why this is (or is not) the case?

(You may also wish to produce a graph for VIEW 2 to help you come up with a good answer. Hint: take a careful look at the three-thread data-point.)

### Answer
1. No
2. It is a load balancing problem. As you can see below, when we use three threads to execute the program, the most computational data points are in the 2nd thread.

   In comparison with **VIEW 2** and other threads, they are not distributed more skewed.

### Result
#### VIEW 1
![mandelbrot-thread-view1.png](https://hackmd.io/_uploads/SkZqUMSQp.png)
![SpeedUp-View1.png](https://hackmd.io/_uploads/r1_Bb1S7p.png)

#### VIEW 2
![mandelbrot-thread-view2.png](https://hackmd.io/_uploads/H135LzSXa.png)
![SpeedUp-View2.png](https://hackmd.io/_uploads/rJPL-JHmT.png)

## Q2
### Question
How do your measurements explain the speedup graph you previously created?

### Answer
The problem is true. Most of computation are in the 2nd thread.

So, that's why the speed up of three-thread is lower than two-thread's.

### Result
#### VIEW 1
![Duration-3-threads-view1.png](https://hackmd.io/_uploads/rJVr5MBQ6.png)

#### VIEW 2
![Duration-3-threads-view2.png](https://hackmd.io/_uploads/H1awcfHm6.png)

## Q3
### Question
In your write-up, describe your approach to parallelization and report the final 4-thread speedup obtained.

### Answer
Compared to the serial implementation, I change the spatial decomposition method from block partition to cyclic partition.

### Result
![Duration-4-threads.png](https://hackmd.io/_uploads/rJpMN7B7a.png)

## Q4
### Question
Now run your improved code with eight threads. Is performance noticeably greater than when running with four threads? Why or why not? (Notice that the workstation server provides 4 cores 4 threads.)

### Answer
No, the workstation can only execute 4 thread in parallel. If we use more threads than the workstation support, it will do context switch among threads leading to increase overhead.

So, the speed up of eight-thread is much lower than four-thread is because of the overhead of context switch.

### Result
![Duration-8-threads.png](https://hackmd.io/_uploads/rJ2jVmrmT.png)
