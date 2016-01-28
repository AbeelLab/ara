Ara
=======

## Building

Like with most scala projects, this tool is shipped with a build.xml file.
Use apache's ant utility to compile this.
The build produces a zip file ara-development.zip, that contains all the necessary jar files.
Personally I don't like that so you should unzip this file.


```shell
$> ant
$> unzip -d build ara-development.zip
```

## Usage

The tool can be run with

```shell
$>java -jar ara.jar <tool> <tool options>
```

More details about the unique-markers tool can be found by:
```shell
$>java -jar build/ara.jar unique-markers
```

## Examples

After building, you can run the examples and tests in the `testing/` directory


