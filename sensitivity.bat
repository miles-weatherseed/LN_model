for /f "eol=; tokens=1,2, 3, 4, 5, 6, 7, 8, 9 delims= " %%a in (rust_params.txt) do (
target\debug\LN_model.exe %%a %%c %%e %%f %%g %%h %%i
target\debug\LN_model.exe %%b %%d %%e %%f %%g %%h %%i
) 