#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

import multiprocessing as mp


class Process(mp.Process):
    counter = 0

    def __init__(self, task_queue, result_queue):
        mp.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue

    def run(self):
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                # None is the Poison pill, means shutdown
                self.task_queue.task_done()
                break
            self.counter += 1
            answer = next_task()
            self.task_queue.task_done()
            self.result_queue.put(answer)
        return self.counter


class WriterProcess(mp.Process):
    counter = 0

    def __init__(self, queue, writer):
        mp.Process.__init__(self)
        self.queue = queue
        self.writer = writer

    def run(self):
        while True:
            result = self.queue.get()
            if result == None:
                # None is the Poison pill, means shutdown
                self.writer.close()
                self.queue.task_done()
                break
            self.counter += 1
            self.writer.write(result)
            self.queue.task_done()
        return self.counter


class ParalellRunner(object):

    def __init__(self, task, reader, writer, task_args=tuple()):
        self.task = task
        self.reader = reader
        self.writer = writer
        self.task_args = task_args
        self.num_reads = 0

    def run(self):
        tasks = mp.JoinableQueue(5000)
        results = mp.JoinableQueue(5000)
        num_procs = mp.cpu_count()

        procs = [Process(tasks, results) for i in xrange(num_procs)]
        for proc in procs:
            proc.start()

        writer_proc = WriterProcess(results, self.writer)
        writer_proc.start()

        for read in self.reader:
            self.num_reads += 1
            tasks.put(self.task(read, *self.task_args))

        # Add a poison pill for each process
        for i in xrange(num_procs):
            tasks.put(None)  # None ends the process

        # Wait for all tasks to run
        tasks.join()
        # Wait for results to be processed
        results.join()
        # Kill Writer subprocess
        results.put(None)
