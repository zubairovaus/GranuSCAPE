# Публикация GranuSCAPE на GitHub

Этот архив подготовлен из `GranuSCAPE(1).zip` (версия `0.1.0rc5`). Код алгоритмов не изменён. Удалены кэши, локальный `test.yaml` и технические отчёты о локальных проверках. Добавлены `RELEASE_CHECKLIST.md` и расширенный `.gitignore`.

## Самый простой способ — через сайт GitHub

1. Войдите на https://github.com/ и нажмите **+ → New repository**.
2. Название: `GranuSCAPE`; описание: `Anatomically guided single-cell morphometry and spatial profiling of fluorescent puncta`.
3. Выберите **Private** для начальной проверки. **Не** включайте автоматическое создание README, .gitignore или лицензии: они уже предусмотрены или требуют решения.
4. Нажмите **Create repository**. На странице нового репозитория выберите **uploading an existing file** (Add file → Upload files).
5. Распакуйте ZIP на компьютере. **Загрузите содержимое папки** `GranuSCAPE_GitHub_ready_v0.1.0rc5`, а не ZIP целиком. GitHub в браузере допускает пакетную загрузку перетаскиванием файлов/папок; если структура папок не сохраняется, воспользуйтесь Git-способом ниже.
6. Укажите комментарий коммита `Initial GranuSCAPE 0.1.0rc5 source release candidate` и подтвердите загрузку.
7. Проверьте `src/granuscape/spatial/nulls.py`, `pyproject.toml`, `README.md`, `tests/` и `src/granuscape/models/` в списке файлов GitHub.
8. Пройдите чек-лист `RELEASE_CHECKLIST.md`, выберите лицензию, проверьте модели `.ilp`, тесты и пример запуска. Затем измените видимость репозитория на **Public**.
9. Создайте GitHub **Release** и Git tag (например, `v0.1.0rc5` для предварительного релиза или `v1.0.0` после завершения проверки), подключите Zenodo и заархивируйте релиз с DOI.

## Более надёжный способ — Git в PowerShell

Установите Git for Windows (https://git-scm.com/downloads/win). Распакуйте ZIP, перейдите в папку проекта и выполните команды, подставив свой логин GitHub:

```powershell
cd "C:\path\to\GranuSCAPE_GitHub_ready_v0.1.0rc5"
git init
git branch -M main
git add .
git status
git commit -m "Initial GranuSCAPE 0.1.0rc5 release candidate"
git remote add origin https://github.com/YOUR_USERNAME/GranuSCAPE.git
git push -u origin main
```

На GitHub предварительно создайте **пустой** приватный репозиторий `GranuSCAPE`. Авторизация при `git push` обычно производится через браузер/Git Credential Manager; пароль аккаунта в командной строке не вводите. Если Git сообщает, что имя и email не заданы, выполните `git config --global user.name "Your Name"` и `git config --global user.email "you@example.org"`, заменив значения на свои.

## Проверка перед открытием

```powershell
python -m pip install -e ".[test]"
python -m pytest -q
```

Для полного запуска необходимы отдельные установленные среды Cellpose и ilastik и адаптация конфигурации `examples/experiment.yaml`. Сам архив не содержит исходные изображения CZI или результаты экспериментов.

**Важное:** отсутствие LICENSE — осознанное временное ограничение, а не рекомендация публиковать без лицензии. Не загружайте в репозиторий личные конфиги, оригинальные изображения или результаты с конфиденциальными сведениями.
