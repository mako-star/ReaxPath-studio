# ReaxPath Studio Web

Local FastAPI workbench for the `reaxnet-lite` streaming analyzer.

## Windows 日常打开

双击项目根目录的 `Open-ReaxFF-Workbench.cmd`。它会检查
`http://127.0.0.1:8767`：服务已运行时直接打开页面，未运行时在后台启动后打开。
重复双击不会重复启动服务。电脑重启后再次双击即可，无需重新分析已保存的任务。
若本地保存着 `acceptance-1000-20261027`，会直接打开这次验收结果。

启动器优先使用项目 `.venv`，其次使用当前机器的
`D:\Users\MECHREVO\miniconda3\python.exe`，最后查找 PATH 中的 Python。
启动日志保存在 `web/.runtime/`。服务在本机后台运行，关闭浏览器后可以重新打开；
更新 Python 源码后需重新启动服务才能加载更新。

## 手动安装与启动

```powershell
cd D:\smiletosvg
python -m pip install -e ".[identifiers]"
cd web
python -m pip install -e .
smiletosvg-web
```

Open `http://127.0.0.1:8000`. For multi-gigabyte trajectories, use the local
path mode so the server reads the source file in place. Chunked upload uses
8 MiB requests, supports resume, and accepts files up to 16 GiB.
