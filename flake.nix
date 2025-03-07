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
          cadiback = nixpkgsFor.${system}.callPackage cadiback-package { };
        in
        {
          inherit
            cadiback
            ;
          default = ganak;
        }
      );
    };
}
