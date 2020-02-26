/*
 * Copyright 2004-2020 Hewlett Packard Enterprise Development LP
 * Other additional copyright holders may be indicated within.
 *
 * The entirety of this work is licensed under the Apache License,
 * Version 2.0 (the "License"); you may not use this file except
 * in compliance with the License.
 *
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


use Path;
use Spawn;
use FileSystem;
use MasonUtils;
use MasonHelp;
use MasonEnv;



proc masonNew(args) throws {
  try! {
    if args.size < 3 {
      masonNewHelp();
      exit();
    } else {
      var vcs = true;
      var show = false;
      var isModuleName = false;
      var moduleName = '';
      var name = '';
      var i = 0;
      for arg in args[2..] {
        i = i + 1;
        if arg == '-h' || arg == '--help' {
          masonNewHelp();
          exit();
        }
        else if arg == '--no-vcs' {
          vcs = false;
        }
        else if arg == '--show' {
          show = true;
        }
        else if arg.startsWith('--moduleName') {
          isModuleName = true;
          if arg.startsWith('--moduleName='){
            var res = arg.split("=");
            moduleName = res[2];
            name = args[2];
            break;
          }
          else {
            moduleName = args[i+2];
            name = args[2];
            break;
          }
        }
        else {
          name = arg;
        }
      }
      if isModuleName then {
        if validatePackageName(moduleName) {
          if isDir(name) {
            throw new owned MasonError("A directory named '" + name + "' already exists");
          }
          InitProject(name, moduleName, isModuleName, vcs, show);
        }
      }
      else {
        if validatePackageName(name) {
          if isDir(name) {
            throw new owned MasonError("A directory named '" + name + "' already exists");
          }
          InitProject(name, moduleName='', isModuleName, vcs, show);
        }
      }
    }
  }
  catch e: MasonError {
    writeln(e.message());
    exit(1);
  }
}

proc validatePackageName(name) throws {
  if name == '' {
    throw new owned MasonError("No package name specified");
  }
  else if !isIdentifier(name) {
    throw new owned MasonError("Bad package name '" + name +
                        "' - only Chapel identifiers are legal package names.\n" +  
                        "Please use mason new <illegal-name> --moduleName <legal-name>");
  }
  else if name.count("$") > 0 {
    throw new owned MasonError("Bad package name '" + name +
                        "' - $ is not allowed in package names");
  }
  else {
    return true;
  }
}

proc InitProject(name, moduleName, isModuleName, vcs, show) throws {
  if vcs {
    gitInit(name, show);
    addGitIgnore(name);
  }
  else {
    mkdir(name);
  }
  // Confirm git init before creating files
  if isDir(name) {
    makeBasicToml(name, path=name);
    makeSrcDir(name);
    if isModuleName then makeModule(name, fileName=moduleName);
    else makeModule(name, fileName=name);
    makeTestDir(name);
    makeExampleDir(name);  
    writeln("Created new library project: " + name);
  }
  else {
    throw new owned MasonError("Failed to create project");
  }
}


proc gitInit(name: string, show: bool) {
  var initialize = "git init -q " + name;
  if show then initialize = "git init " + name;
  runCommand(initialize);
}

proc addGitIgnore(name: string) {
  var toIgnore = "target/\nMason.lock\n";
  var gitIgnore = open(name+"/.gitignore", iomode.cw);
  var GIwriter = gitIgnore.writer();
  GIwriter.write(toIgnore);
  GIwriter.close();
}


proc makeBasicToml(name: string, path: string) {
  const baseToml = '[brick]\n' +
                     'name = "' + name + '"\n' +
                     'version = "0.1.0"\n' +
                     'chplVersion = "' + getChapelVersionStr() + '"\n' +
                     '\n' +
                     '[dependencies]' +
                     '\n';
  var tomlFile = open(path+"/Mason.toml", iomode.cw);
  var tomlWriter = tomlFile.writer();
  tomlWriter.write(baseToml);
  tomlWriter.close();
}

proc makeSrcDir(path:string) {
  mkdir(path + "/src");
}

proc makeModule(path:string, fileName:string) {
  const libTemplate = '/* Documentation for ' + fileName +
  ' */\nmodule '+ fileName + ' {\n  writeln("New library: '+ fileName +'");\n}';
  var lib = open(path+'/src/'+fileName+'.chpl', iomode.cw);
  var libWriter = lib.writer();
  libWriter.write(libTemplate + '\n');
  libWriter.close();
}

proc makeTestDir(path:string) {
  mkdir(path + "/test");
}

proc makeExampleDir(path:string) {
  mkdir(path + "/example");
}
