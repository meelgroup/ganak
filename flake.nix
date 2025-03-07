{
  description = "A Scalable Probabilistic Exact Model Counter";
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixpkgs-unstable";
  };

  outputs =
    {
      self,
      nixpkgs,
    }:
    let
      inherit (nixpkgs) lib;
      systems = lib.intersectLists lib.systems.flakeExposed lib.platforms.linux;
      forAllSystems = lib.genAttrs systems;
      nixpkgsFor = forAllSystems (system: nixpkgs.legacyPackages.${system});

      cadical-package =
        {
          stdenv,
          fetchFromGitHub,
          lsd,
        }:
        stdenv.mkDerivation {
          name = "cadical";
          src = fetchFromGitHub {
            owner = "meelgroup";
            repo = "cadical";
            rev = "45850b35836122622a983e637251299cc16f3161";
            name = "cadical";
            hash = "sha256-ugAudDgw91DeR90zQR5RlCKoLv/hDdv03oa1v3lG1nY=";
          };

          configurePhase = ''./configure --competition'';

          installPhase = ''
            mkdir -p $out/lib
            cp -r . $out
            cp build/libcadical.a $out/lib
            mkdir -p $out/include
            cp src/*.hpp $out/include
          '';
        };
      cadiback-package =
        {
          stdenv,
          fetchFromGitHub,
          git,
          cadical,
        }:
        stdenv.mkDerivation {
          name = "cadiback";
          src = fetchFromGitHub {
            owner = "meelgroup";
            repo = "cadiback";
            rev = "860b1df7f7b5e966b423163d89ec5a972e09ae67";
            name = "cadiback";
            hash = "sha256-sKKIOWQjW2c5Duz/jf00ggwXMlQ9KL7mVL7val47Cng=";
            leaveDotGit = true;
          };

          nativeBuildInputs = [ git ];
          buildInputs = [ cadical ];

          patchPhase = ''
            substituteInPlace makefile.in \
              --replace-fail "/usr/" "$out/" \
              --replace-fail "../cadical" "${cadical}"
          '';
          configurePhase = ''
            export CADICAL=${cadical} 
            ./configure
          '';

          preInstall = ''
            mkdir -p $out/lib
            mkdir -p $out/include
          '';
        };
      ensmallen-package =
        {
          stdenv,
          cmake,
          fetchzip,
          armadillo,
        }:
        stdenv.mkDerivation {
          name = "ensmallen";
          src = fetchzip {
            url = "https://ensmallen.org/files/ensmallen-2.21.1.tar.gz";
            hash = "sha256-6LZooaR0rmqWgEm0AxmWoVPuIahjOfwSFu5cssc7LA8=";
          };
          nativeBuildInputs = [ cmake ];
          buildInputs = [ armadillo ];
        };
      mlpack-package =
        {
          stdenv,
          fetchFromGitHub,
          lsd,
          cmake,
          armadillo,
          ensmallen,
          cereal,
          pkg-config,
        }:
        stdenv.mkDerivation {
          name = "mlpack";
          src = fetchFromGitHub {
            "owner" = "mlpack";
            "repo" = "mlpack";
            "rev" = "4.4.0";
            "hash" = "sha256-EPz8qPTUAldS+k5/qkZf8EKXKjnxElfJxlTEMLPhTQE=";
          };
          nativeBuildInputs = [
            pkg-config
            cmake
            armadillo
          ];
          buildInputs = [
            ensmallen
            cereal
          ];
        };
      sbva-package =
        {
          stdenv,
          eigen,
          fetchFromGitHub,
          cmake,
          autoPatchelfHook,
        }:
        stdenv.mkDerivation {
          name = "sbva";
          src = fetchFromGitHub {
            owner = "meelgroup";
            repo = "SBVA";
            rev = "5912435affe8c77ecf364093cea29e0fc5c1b5cb";
            hash = "sha256-BoR14FBH3eCPYio6l6d+oQp3/hu4U7t1STb9NgSWJ2M=";
          };
          nativeBuildInputs = [
            cmake
            autoPatchelfHook
          ];
          buildInputs = [ eigen ];
        };
      breakid-package =
        {
          stdenv,
          cmake,
          autoPatchelfHook,
          fetchFromGitHub,
        }:
        stdenv.mkDerivation {
          name = "breakid";
          src = fetchFromGitHub {
            owner = "meelgroup";
            repo = "breakid";
            rev = "f90eff992342024640f9d539c1e68ef924b22994";
            hash = "sha256-hxVvfOfPQxZ6RUYv6BWWohEbvoPfjWIa+UNtkOuvS1Q=";
          };
          nativeBuildInputs = [
            cmake
            autoPatchelfHook
          ];
          buildInputs = [ ];
        };
    in
    {
      # checks

      devShells = forAllSystems (
        system:
        let
          pkgs = nixpkgsFor.${system};
        in
        {
          default = pkgs.mkShell {
            packages = [
              pkgs.cmake
              pkgs.pkg-config
              pkgs.gmp
            ];
          };
        }
      );

      packages = forAllSystems (
        system:
        let
          ganak-package =
            {
              stdenv,
              cmake,
              pkg-config,
              gmp,
              autoPatchelfHook,
            }:
            stdenv.mkDerivation {
              name = "ganak";
              nativeBuildInputs = [
                cmake
                pkg-config
                autoPatchelfHook
              ];
              buildInputs = [ gmp ];
              src = ./.;
              installPhase = ''
                mkdir -p $out/bin $out/lib
                mv ./ganak $out/bin/
                mv ./src/libganak.so.* \
                   ./src/clhash/libclhash.so \
                   ./src/component_types/libcomponent_types.so \
                   $out/lib
              '';
            };
          ganak = nixpkgsFor.${system}.callPackage ganak-package { };
          cadical = nixpkgsFor.${system}.callPackage cadical-package { };
          cadiback = nixpkgsFor.${system}.callPackage cadiback-package { inherit cadical; };
          sbva = nixpkgsFor.${system}.callPackage sbva-package { };
          breakid = nixpkgsFor.${system}.callPackage breakid-package { };
          ensmallen = nixpkgsFor.${system}.callPackage ensmallen-package { };
          mlpack = nixpkgsFor.${system}.callPackage mlpack-package { inherit ensmallen; };
        in
        {
          inherit
            cadical
            sbva
            breakid
            cadiback
            ensmallen
            mlpack
            ;
          default = ganak;
        }
      );
    };
}
