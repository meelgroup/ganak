// Runs one Emscripten CLI tool (ganak.js or bdnnf.js) per worker: write the
// input files, call main, return stdout/stderr, the exit status and the output files.

interface RunRequest {
  version: string;
  tool: string;
  args: string[];
  files: Record<string, string | Uint8Array>;
  outFiles: string[];
}

self.onmessage = (e: MessageEvent<RunRequest>) => {
  const req = e.data;
  const out: string[] = [];
  const err: string[] = [];
  const mod: any = {
    noInitialRun: true,
    // the build's version stamp in every URL: a new build is never served from cache
    locateFile: (p: string) => p + "?v=" + encodeURIComponent(req.version),
    print: (t: string) => out.push(t),
    printErr: (t: string) => {
      if (!t.includes("__syscall_getrusage")) err.push(t);
    },
    onRuntimeInitialized: () => {
      for (const [p, d] of Object.entries(req.files)) mod.FS.writeFile(p, d);
      let status = 0;
      try {
        const r = mod.callMain(req.args);
        if (typeof r === "number") status = r;
      } catch (x: any) {
        // exit() inside the tool throws ExitStatus
        if (x && typeof x.status === "number") status = x.status;
        else {
          err.push(String(x));
          status = -1;
        }
      }
      const files: Record<string, Uint8Array> = {};
      for (const p of req.outFiles) {
        try {
          files[p] = mod.FS.readFile(p);
        } catch {
          // not produced
        }
      }
      postMessage({ status, stdout: out.join("\n"), stderr: err.join("\n"), files });
    },
  };
  (self as any).Module = mod;
  importScripts(req.tool + ".js?v=" + encodeURIComponent(req.version));
};
