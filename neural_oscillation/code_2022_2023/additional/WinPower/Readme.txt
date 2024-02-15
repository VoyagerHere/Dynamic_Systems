WinPower - Mex function to set the power status on Windows computers
Power off, shut down, logoff, sleep, reboot, reboot and restart Matlab, disable/enable sleeping, lock
WinPower(Command, Force)
INPUT:
Command: String, not case-sensitive.
'poweroff': Switch power off.
'reboot': Reboot the machine.
'logoff': Logoff the current user.
'shutdown': Shut down the machine to a state, which allows the user to
switch off the power securely by hand.
'sleep': Let the machine fall asleep.
'sleep', 'off': Disable the sleep timers, 'on' enables them.
'hibernate': Write memory to disk and fall into deep sleep.
'lock': Lock the machine, password is required for wake-up.
'rebootmatlab': Reboot the machine, restart Matlab when the user is
logged in again.
'monitor': Set monitor status without stopping the processing.
2nd input: 'off' (default), 'on', 'standby'.
Moving the mouse etc. enables the monitor automatically.
'battery': Reply the battery related parameters for laptops.
Force: Optional argument to force the action:
'force': Close waiting applications. Dangerous.
'forceifhung': Close waiting and crashed applications. Dangerous.

EXAMPLES:
1. Try a poweroff, do not kill waiting applications:
WinPower('poweroff');
2. Force a poweroff (unsaved documents are lost!):
WinPower('poweroff', 'forceifhung');
3. Do not let the computer fall asleep during a long computation:
WinPower('Sleep', 'off'); Long_Calculation(); WinPower('Sleep', 'on');
4. Reboot the machine and after the user is logged in Matlab is started with
the bench() function:
WinPower('RebootMatlab', 'bench(2)');
5. Get the battery parameters:
Status = WinPower('BatteryStatus')
6. Switch off the monitor for 5 seconds:
WinPower('Monitor', 'off'); pause(5); WinPower('Monitor', 'on');

WinPower runs under Windows only. It can be compiled with MSVC2008/2010, but LCC shipped with Matlab 32bit fails. In case of troubles or for Matlab 6.5 download pre-compiled files:
http://www.n-simon.de/mex