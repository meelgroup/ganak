// Headless (Node.js) test of the in-browser d-DNNF pipeline. Runs the same chain
// the page does, without a browser:
//   ganak --compile (raw BDNNF) -> bdnnf cleanup (clean BDNNF) -> bdnnf dot (DOT)
// using the same Emscripten modules shipped here. A smoke test after rebuilding wasm.
//
// These are non-MODULARIZE Emscripten CLI modules expecting a global `Module`. The
// browser gives each its own Web Worker; here we run each in a fresh function scope
// with `Module` in scope (a plain `require` would let its own `var Module` shadow ours).
//
// Usage:
//   node chain.js [input.cnf] [out.dot]
//   node chain.js                      # built-in pwmc example, prints DOT

const fs = require('fs');
const path = require('path');

const HERE = __dirname;

// Default example matches the page's textarea (a projected weighted CNF).
const DEFAULT_CNF = `c t pwmc
p cnf 3 2
c p weight 1 4/10 0
c p weight 2 6/10 0
c p weight 3 8/10 0
1 2 0
-3 -1 -2 0
c p show 1 2 3 0
`;

// Run one Emscripten CLI tool: write inFiles into its in-memory FS, call main with
// args, return outFile's contents (or null). Mirrors the worker in index.html.
function runTool(jsFile, args, inFiles, outFile) {
  return new Promise((resolve, reject) => {
    const code = fs.readFileSync(path.join(HERE, jsFile), 'utf8');
    const Module = {
      noInitialRun: true,
      arguments: [],
      print: t => process.stderr.write('[' + jsFile + '] ' + t + '\n'),
      printErr: t => {
        if (!t.includes('__syscall_getrusage')) process.stderr.write('[' + jsFile + '] ' + t + '\n');
      },
      onRuntimeInitialized: function () {
        try {
          for (const p in inFiles) Module.FS.writeFile(p, inFiles[p]);
          Module.callMain(args);
          resolve(!outFile ? null : outFile.endsWith('.bdnnf') ? Module.FS.readFile(outFile)
                                   : Module.FS.readFile(outFile, { encoding: 'utf8' }));
        } catch (e) { reject(e); }
      }
    };
    // Keep the tool from auto-running main off our argv.
    process.argv = [process.argv[0]];
    new Function('Module', 'require', '__dirname', '__filename', code)(
      Module, require, HERE, path.join(HERE, jsFile));
  });
}

(async () => {
  const inPath = process.argv[2];
  const outPath = process.argv[3];
  const cnf = inPath ? fs.readFileSync(inPath, 'utf8') : DEFAULT_CNF;

  // 1. ganak --compile: the raw binary BDNNF circuit (faithful by default).
  //    --mode 1 (weighted rational) handles both weighted and plain CNFs.
  const raw = await runTool('ganak.js',
    ['--verb', '1', '--mode', '1', '--compile', '/out.bdnnf', '/input.cnf'],
    { '/input.cnf': cnf }, '/out.bdnnf');
  process.stderr.write('\n=== RAW BDNNF (' + raw.length + ' bytes) ===\n');

  // 2. bdnnf cleanup: flatten, drop dead nodes.
  const cleaned = await runTool('bdnnf.js',
    ['cleanup', '/in.bdnnf', '/out.bdnnf'], { '/in.bdnnf': raw }, '/out.bdnnf');

  // 3. bdnnf dot: render the cleaned circuit as Graphviz DOT.
  const dot = await runTool('bdnnf.js',
    ['dot', '/in.bdnnf', '/out.dot'], { '/in.bdnnf': cleaned }, '/out.dot');

  if (outPath) {
    fs.writeFileSync(outPath, dot);
    process.stderr.write('\nwrote DOT to ' + outPath + '\n');
  } else {
    process.stdout.write(dot);
  }
  process.exit(0);
})().catch(e => { process.stderr.write('FAIL: ' + e + '\n'); process.exit(1); });
