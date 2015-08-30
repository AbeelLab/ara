package macaw2

abstract class Marker {

}

object Marker {
  def unapply(s: String): Option[(String, Int, Boolean)] = {
    val line = s.split("\t")
    val mutInfo = line(0)
    val count = line(1).toInt
    val presence = if (line(3) == "present") true else false
    Some(mutInfo, count, presence)
  }
}

case class DrugMarker(val mutInfo: String, val count: Int, val isPresent: Boolean) extends Marker {
  override def toString(): String = mutInfo + "\t" + count + "\t" + (if (isPresent) "present" else "absent")
  val arr = mutInfo.split("_")
  val coordinate = arr(0).drop(1).dropRight(1).toInt
  val markerType = arr(1)
  val drugs = arr(2)
}

case class ClusterMarker(val mutInfo: String, val count: Int, val isPresent: Boolean) extends Marker {
  val lineage = mutInfo.split("_")(2)
  override def toString(): String = mutInfo + "\t" + count + "\t" + (if (isPresent) "present" else "absent") + "\t" + lineage
}

