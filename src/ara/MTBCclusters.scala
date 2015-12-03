package ara

trait MTBCclusters {

  val cFile = scala.io.Source.fromInputStream(AraUtilities.getClass().getResourceAsStream("/hierclusters_global_numbers.txt")).getLines().filterNot(f => f.startsWith("#") || f.trim.size == 0).toList
  
  val root = Map("MTBC" -> Map("size" -> 5993, "markers" -> 0))
  val cNumbers = cFile.map { line =>
    val arr = line.split("\t")
    (arr(0), Map("size" -> arr(1).toInt, "markers" -> (arr(8).toInt + arr(9).toInt)))
  }.toMap ++ root

  val mtbcClusters = cNumbers.keysIterator.toList

}