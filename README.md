# Practice in MAGE
praktika ;)

Сделать обработчик изображений с помощью перцептивного хэша.

Использованный библиотеки: CImg, (и заголовочный файл pHash)

Нужно было сделать:

1)Чтение какой-то директории на наличие изображений, если имеются, то переходим к шагу 2, если же нет, то exit(1);

2)Чтение файлов и последующая обработка в соотвествии с алгоритмом перцептивного хэша:

  a) ресайс изображения(32x32)
  
  b)  обесцвечивание
  
  c) перевод битов в матрицу 8x8
  
  d) подсчет хэша
  
3)После выполнения шага 2, создание xml файла с именем файла изображения, куда записали бы: оригинальный размер файла, имя, хэш

4)Обработка таким способом всех изображений

5)Последующее их сравнение и разбитие в папки по индексам(схожие)

6)Написать интерефейс(не знаю как ;D)


Что сделано:

Первый, второй пункты.



