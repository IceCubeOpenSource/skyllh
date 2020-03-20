# -*- coding: utf-8 -*-

import numpy as np
import time

from skyllh.core import display
from skyllh.core.py import range

"""The timing module provides code execution timing functionalities. The
TimeLord class keeps track of execution times of specific code segments,
called "tasks". The TaskTimer class can be used within a `with`
statement to time the execution of the code within the `with` block.
"""


class TaskRecord(object):
    def __init__(self, name, tstart, tend):
        """Creates a new TaskRecord instance.

        Parameters
        ----------
        name : str
            The name of the task.
        tstart : float | 1d ndarray of float
            The start time(s) of the task in seconds.
        tend : float | 1d ndarray of float
            The end time(s) of the task in seconds.
        """
        self.name = name

        tstart = np.atleast_1d(tstart)
        tend = np.atleast_1d(tend)

        if(len(tstart) != len(tend)):
            raise ValueError('The number of start and end time stamps must '
                             'be equal!')

        # Create a (2,Niter)-shaped 2D ndarray holding the start and end time
        # stamps of the task executions. This array must be sorted by the start
        # time stamps.
        self._tstart_tend_arr = np.sort(np.vstack((tstart, tend)), axis=1)

    @property
    def tstart(self):
        """(read-only) The time stamps the execution of this task started.
        """
        return self._tstart_tend_arr[0, :]

    @property
    def tend(self):
        """(read-only) The time stamps the execution of this task was stopped.
        """
        return self._tstart_tend_arr[1, :]

    @property
    def duration(self):
        """(read-only) The total duration (without time overlap) the task was
        executed.
        """
        arr = self._tstart_tend_arr

        d = arr[1, 0] - arr[0, 0]
        last_tend = arr[1, 0]
        n = self.niter
        for idx in range(1, n):
            tstart = arr[0, idx]
            tend = arr[1, idx]
            if(tend <= last_tend):
                continue
            if(tstart <= last_tend and tend > last_tend):
                d += tend - last_tend
            elif(tstart >= last_tend):
                d += tend - tstart
            last_tend = tend

        return d

    @property
    def niter(self):
        """(read-only) The number of times this task was executed.
        """
        return self._tstart_tend_arr.shape[1]

    def join(self, tr):
        """Joins this TaskRecord with the given TaskRecord instance.

        Parameters
        ----------
        tr : TaskRecord
            The instance of TaskRecord that should be joined with this
            TaskRecord instance.
        """
        self._tstart_tend_arr = np.sort(
            np.append(
                self._tstart_tend_arr, np.vstack((tr.tstart, tr.tend)),
                axis=1),
            axis=1)


class TimeLord(object):
    def __init__(self):
        self._task_records = []
        self._task_records_name_idx_map = {}

    @property
    def task_name_list(self):
        """(read-only) The list of task names.
        """
        return list(self._task_records_name_idx_map.keys())

    def add_task_record(self, tr):
        """Adds a given task record to the internal list of task records.
        """
        tname = tr.name

        if(self.has_task_record(tname)):
            # The TaskRecord already exists. Update the task record.
            self_tr = self.get_task_record(tname)
            self_tr.join(tr)
            return

        self._task_records.append(tr)
        self._task_records_name_idx_map[tr.name] = len(self._task_records)-1

    def get_task_record(self, name):
        """Retrieves a task record of the given name.

        Parameters
        ----------
        name : str
            The name of the task record.

        Returns
        -------
        task_record : TaskRecord
            The instance of TaskRecord with the requested name.
        """
        return self._task_records[self._task_records_name_idx_map[name]]

    def has_task_record(self, name):
        """Checks if this TimeLord instance has a task record of the given name.

        Parameters
        ----------
        name : str
            The name of the task record.

        Returns
        -------
        check : bool
            ``True`` if this TimeLord instance has a task record of the given
            name, and ``False`` otherwise.
        """
        return name in self._task_records_name_idx_map

    def join(self, tl):
        """Joins a given TimeLord instance with this TimeLord instance. Tasks
        of the same name will be updated and new tasks will be added.

        Parameters
        ----------
        tl : TimeLord instance
            The instance of TimeLord whos tasks should be joined with the tasks
            of this TimeLord instance.
        """
        for tname in tl.task_name_list:
            other_tr = tl.get_task_record(tname)
            if(self.has_task_record(tname)):
                # Update the task record.
                tr = self.get_task_record(tname)
                tr.join(other_tr)
            else:
                # Add a new task record.
                self.add_task_record(other_tr)

    def task_timer(self, name):
        """Creates TaskTimer instance for the given task name.
        """
        return TaskTimer(self, name)

    def __str__(self):
        """Generates a pretty string for this time lord.
        """
        s = 'Executed tasks:'
        task_name_list = self.task_name_list
        task_name_len_list = [len(task_name) for task_name in task_name_list]
        max_task_name_len = np.minimum(
            np.max(task_name_len_list), display.PAGE_WIDTH-25)

        n_tasks = len(task_name_list)
        for i in range(n_tasks):
            tr = self._task_records[i]
            task_name = tr.name[0:max_task_name_len]
            t = tr.duration / tr.niter
            line = '\n[{task_name:'+str(max_task_name_len)+'s}] {t:7.{p}{c}} '\
                'sec/iter ({niter:d})'
            s += line.format(
                task_name=task_name,
                t=t,
                p=1 if t > 1e3 or t < 1e-3 else 3,
                c='e' if t > 1e3 or t < 1e-3 else 'f',
                niter=tr.niter)

        return s


class TaskTimer(object):
    def __init__(self, time_lord, name):
        """
        Parameters
        ----------
        time_lord : instance of TimeLord
            The TimeLord instance that keeps track of the recorded tasks.
        name : str
            The name of the task.
        """
        self.time_lord = time_lord
        self.name = name

        self._start = None
        self._end = None

    @property
    def time_lord(self):
        """The TimeLord instance that keeps track of the recorded tasks. This
        can be None, which means that the task should not get recorded.
        """
        return self._time_lord

    @time_lord.setter
    def time_lord(self, lord):
        if(lord is not None):
            if(not isinstance(lord, TimeLord)):
                raise TypeError('The time_lord property must be None or an '
                                'instance of TimeLord!')
        self._time_lord = lord

    @property
    def name(self):
        """The name if the task.
        """
        return self._name

    @name.setter
    def name(self, name):
        if(not isinstance(name, str)):
            raise TypeError('The name property must be an instance of str!')
        self._name = name

    @property
    def duration(self):
        """The duration in seconds the task was executed.
        """
        return (self._end - self._start)

    def __enter__(self):
        """This gets executed when entering the `with` block.
        """
        self._start = time.time()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """This gets executed when exiting the `with` block.
        """
        self._end = time.time()

        if(self._time_lord is None):
            return

        self._time_lord.add_task_record(TaskRecord(
            self._name, self._start, self._end))
