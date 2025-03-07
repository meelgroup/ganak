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

      cadiback-package =
        {
          stdenv,
          fetchFromGitHub,
          git,
        }:
        stdenv.mkDerivation {
          name = "cadiback";
          srcs = [
            (fetchFromGitHub {
              owner = "meelgroup";
              repo = "cadiback";
              rev = "860b1df7f7b5e966b423163d89ec5a972e09ae67";
              name = "cadiback";
              hash = "sha256-sKKIOWQjW2c5Duz/jf00ggwXMlQ9KL7mVL7val47Cng=";
              leaveDotGit = true;
            })
            (fetchFromGitHub {
              owner = "meelgroup";
              repo = "cadical";
              rev = "45850b35836122622a983e637251299cc16f3161";
              name = "cadical";
              hash = "sha256-ugAudDgw91DeR90zQR5RlCKoLv/hDdv03oa1v3lG1nY=";
            })
          ];
          sourceRoot = "cadiback";

          nativeBuildInputs = [ git ];

          patchPhase = ''
            substituteInPlace makefile.in \
              --replace-fail "/usr/" "$out/"
          '';
          # I don't like building cadical in this phase, but cadiback requires the build outputs during configuring.
          configurePhase = ''
            chmod -R u+w ../cadical
            ln -s ../cadical
            mkdir ../cadical/build
            cd ../cadical
            ./configure --competition && make
            cd ../cadiback
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
          sbva = nixpkgsFor.${system}.callPackage sbva-package { };
          cadiback = nixpkgsFor.${system}.callPackage cadiback-package { };
          breakid = nixpkgsFor.${system}.callPackage breakid-package { };
          ensmallen = nixpkgsFor.${system}.callPackage ensmallen-package { };
        in
        {
          inherit
            sbva
            breakid
            cadiback
            ensmallen
            ;
          default = ganak;
        }
      );
    };
}
