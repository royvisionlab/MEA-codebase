@echo off

for /d /R %%D in (..\data\20220412C\*) do echo %%~fD

for /d /R %%D in (..\data\20220412C\*) do python convert_litke_to_kilosort.py %%~fD C:\Users\manoo\Documents\SpikeSorting\kilosort %%D -w -v 3





