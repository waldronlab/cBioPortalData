The cBioPortal service uses OpenAPI 3.0, but our software (currently)
uses rapiclient, with a restriction to Swagger 2.0. The file in this
directory was created using the [LucyBot][] api-spec-converter from
the command line.

## Setup

```
apt install nodejs
apt install npm
```

## Conversion software installation

```
npm install -g @apiture/openapi-down-convert
npm install -g api-spec-converter
```

## API Protocol Version Conversion

```
## download the OpenAPI 3.1 specification
wget -O openapi.json https://www.cbioportal.org/api/v3/api-docs/internal
## convert 3.1 to 3.0
openapi-down-convert --input openapi.json --output api3.0.json
## convert 3.0 to 2.0
api-spec-converter -f openapi_3 -t swagger_2 \
  api3.0.json > \
  api.json
```

[LucyBot]: https://github.com/LucyBot-Inc/api-spec-converter
